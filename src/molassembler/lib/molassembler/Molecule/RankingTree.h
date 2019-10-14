/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief IUPAC-like ranking of substituents at atoms
 *
 * Centerpoint of library ranking algorithm. Implements the RankingTree class,
 * which can be instantiated on any atomic index in a Molecule, which
 * splits the Molecule into an acyclic tree. That tree can then be used to rank
 * that atom's direct substituents according to IUPAC-like sequence rules.
 *
 * @todo
 * - Rewrite the shebang with Posets at nodes for transferability and
 *   persistence of ranking
 * - Consider transition to more unordered containers
 * - Pseudo-asymmetry considerations
 *   - Is the propagation of pseudoasymmetry REALLY necessary? It makes sense
 *     to permute a stereopermutator designated as pseudo-asymmetric. I don't
 *     know if creating the exact diastereomer of a molecule will ever be
 *     needed as an operation (if so, then yes, pseudoasymmetry pop is needed).
 *     Creating ALL stereomers of a molecule can be done without
 *     pseudoasymmetry.  All other points below are only relevant in case
 *     pseudo-asymmetry needs to be included.
 *   - Stereopermutator interface change to support pseudo-asymmetry tag (?)
 *   - Ranking function interface change to propagate pseudo-asymmetry result
 * - When in _auxiliaryApplySequenceRules, and BFS is omnidirectional, isn't
 *   the up-edge always top-ranked? It may not be necessary to include it in
 *   the comparisonSets and as a result, also not complicate BFS handling.
 */

#ifndef INCLUDE_MOLASSEMBLER_RANKING_HIERARCHICAL_TREE_H
#define INCLUDE_MOLASSEMBLER_RANKING_HIERARCHICAL_TREE_H

#include "boost/variant.hpp"

#include "molassembler/Detail/BuildTypeSwitch.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/Log.h"
#include "molassembler/Molecule.h"
#include "molassembler/Containers/OrderDiscoveryHelper.h"

#include "temple/Adaptors/AllPairs.h"

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

namespace Scine {

namespace molassembler {

/*!
 * @brief Central class for unified IUPAC-like ranking of organic and inorganic
 *   structures.
 *
 * The general procedure is that any cycles are broken down in a BFS
 * iteration through the molecular graph at the atom vertex of instantiation.
 * Afterwards, all branches are expanded in a DFS-like manner to ensure that
 * cycle closures are properly reported as duplicate atoms. In the rank() member
 * function, the actual ranking as demanded by the IUPAC sequence rules (adapted
 * for the mixed inorganic/organic molecular graphs) is performed.
 */
class RankingTree {
private:
  /*! For properly naming all the log files emitted in case the appropriate log
   * particular is set.
   */
  static unsigned _debugMessageCounter;

  /*! Modifies the graphviz molgraph into a digraph so that it can be combined
   * with all the other digraphs using gvpack
   */
  static std::string _adaptMolGraph(std::string molGraph);

  //! Writes graphviz log files from graph strings
  static void _writeGraphvizFiles(
    const std::vector<std::string>& graphvizStrings
  );

public:
//!@name Member types
//!@{
  /*! Option deciding in what manner the ranking tree is expanded
   *
   * The optimized expansion only expands positions that are needed for ranking
   * while applying sequence rule 1 immediately. Full expansion expands the
   * entire tree indiscriminately of whether the information will ever be needed
   * prior to applying any sequence rules.
   */
  enum class ExpansionOption {
    //! Expand the tree only as needed
    OnlyRequiredBranches,
    //! Exhaustively expand all vertices until completion before ranking
    Full
  };

  //! Data class that sets which supplementary data is stored for a tree vertex
  struct VertexData {
    //! The original molecule index this vertex represents
    AtomIndex molIndex;
    //! Whether this vertex in the tree represents a duplicate
    bool isDuplicate;

    //! A vertex-central stereopermutator may be instantiated here or may not
    boost::optional<AtomStereopermutator> stereopermutatorOption;

    // This is the place for eventual atomic number deviations
    // boost::optional<double> deviantAtomicNumber;
  };

  //! Data class that sets which supplementary data is stored for a tree edge
  struct EdgeData {
    boost::optional<BondStereopermutator> stereopermutatorOption;
  };

  //! The BGL Graph type used to store the tree
  using BGLType = boost::adjacency_list<
    /* OutEdgeListS = Type of Container for edges of a vertex
     * Options: vector, list, slist, set, multiset, unordered_set
     * Choice: setS
     */
    boost::setS,
    /* VertexListS = Type of Container for vertices
     * Options: vector, list, slist, set, multiset, unordered_set
     * Choice: vecS, removing vertices does not occur
     */
    boost::vecS,
    /* DirectedS = Is the graph directed or not?
     * Choice: Bidirectional. Although a tree is merely directional, and not
     *   bi-directional, we nevertheless also need fast access to in_edges, so
     *   we choose bidirectional
     */
    boost::bidirectionalS,
    // VertexProperty = What information is stored about vertices?
    VertexData,
    // EdgeProperty = What information is stored about edges?
    EdgeData
  >;

  // IUPAC Sequence rule one tree vertex comparator
  class SequenceRuleOneVertexComparator;

  // IUPAC Sequence rule three tree edge comparator
  class SequenceRuleThreeEdgeComparator;

  // IUPAC Sequence rule four tree vertex and edge mixed comparator
  class SequenceRuleFourVariantComparator;

  // IUPAC Sequence rule five tree vertex and edge mixed comparator
  class SequenceRuleFiveVariantComparator;

  // Returns whether the vertex or edge has an instantiated Stereopermutator
  class VariantHasInstantiatedStereopermutator;

  // Returns the mixed depth (see _mixedDepth functions) of the vertex or edge
  class VariantDepth;

  // Returns the source node of a vertex (identity) or an edge (source node)
  class VariantSourceNode;

  // Returns a string representation of the stereopermutator on vertex or edge
  class VariantStereopermutatorStringRepresentation;

  /*
   * Returns whether two variants' stereopermutators form a like pair or not
   * according to the sub-rules in IUPAC sequence rule four.
   */
  class VariantLikePair;

  //! Helper class to write a graphviz representation of the generated tree
  class GraphvizWriter;

  //! Data class to store junction vertex and paths from the source vertices
  struct JunctionInfo;
//!@}

private:
  //! Type of tree vertex accessor
  using TreeVertexIndex = BGLType::vertex_descriptor;
  //! Type of tree edge accessor
  using TreeEdgeIndex = BGLType::edge_descriptor;
  //! Variant type of both
  using VariantType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

  //! A readbility-improving constexpr replacement for the root index 0 in code
  static constexpr TreeVertexIndex rootIndex = 0;

/* State */
  //! The BGL Graph representing the acyclic tree
  BGLType _tree;

  //! The helper instance for discovering the ordering of the to-rank branches
  OrderDiscoveryHelper<TreeVertexIndex> _branchOrderingHelper;

  // Closures
  const OuterGraph& _graphRef;
  const Cycles& _cyclesRef;
  const StereopermutatorList& _stereopermutatorsRef;
  const std::string _adaptedMolGraphviz;

/* Minor helper classes and functions */
  //! Returns the parent of a node. Fails if called on the root!
  TreeVertexIndex _parent(const TreeVertexIndex& index) const;

  //! Returns an unordered list of all adjacent vertices (in- and out-adjacents)
  std::vector<TreeVertexIndex> _adjacents(TreeVertexIndex index) const;

  //! Counts the number of terminal hydrogens on out-edges of the specified index
  unsigned _adjacentTerminalHydrogens(const TreeVertexIndex& index) const;

  //! Checks whether a tree index is the result of a bond split
  bool _isBondSplitDuplicateVertex(const TreeVertexIndex& index) const;

  //! Checks whether a tree index is the result of a cycle closure
  bool _isCycleClosureDuplicateVertex(const TreeVertexIndex& index) const;

  /*! Determines which tree indices are to be ranked
   *
   * Determines which tree indices are to be ranked based on the central index
   * and a list of index excludes
   */
  std::vector<TreeVertexIndex> _auxiliaryAdjacentsToRank(
    TreeVertexIndex sourceIndex,
    const std::vector<TreeVertexIndex>& excludeIndices
  ) const;

  //! Returns whether any of the comparison multisets has an entry
  template<typename ComparisonSets>
  bool _notEmpty(const ComparisonSets& comparisonSets) const {
    return temple::all_of(
      comparisonSets,
      [](const auto& mapPair) -> bool {
        return !mapPair.second.empty();
      }
    );
  }

  /*!
   * Returns the node degree, excluding edges to or from nodes marked as
   * duplicate
   */
  unsigned _nonDuplicateDegree(const TreeVertexIndex& index) const;

  /*! Adds the required duplicate atoms for bonds of high bond orders
   *
   * This adds the appropriate amount of duplicate atoms to the source and
   * target vertices according to the bond order (if integral). Returns
   * the new duplicate tree vertex indices of duplicate nodes added to the
   * source vertex only, not those added to the target vertex.
   */
  std::vector<TreeVertexIndex> _addBondOrderDuplicates(
    const TreeVertexIndex& treeSource,
    const TreeVertexIndex& treeTarget
  );

  //! Returns all tree indices in the branch from the specified index up to root
  std::unordered_set<TreeVertexIndex> _treeIndicesInBranch(TreeVertexIndex index) const;

  //! Returns all mol indices in the branch from the specified index up to root
  std::unordered_set<AtomIndex> _molIndicesInBranch(TreeVertexIndex index) const;

  //! Returns a depth measure of a duplicate vertex for sequence rule 1
  unsigned _duplicateDepth(TreeVertexIndex index) const;

  //! Returns the depth of a node in the tree
  unsigned _depthOfNode(TreeVertexIndex index) const;

  //! Returns a mixed depth measure for ranking both vertices and edges
  unsigned _mixedDepth(const TreeVertexIndex& vertexIndex) const;

  //! Returns a mixed depth measure for ranking both vertices and edges
  unsigned _mixedDepth(const TreeEdgeIndex& edgeIndex) const;

  /*!
   * Returns the deepest vertex in the tree in whose child branches both a and
   * b are located
   */
  JunctionInfo _junction(const TreeVertexIndex& a, const TreeVertexIndex& b) const;

  //! Returns whether a molecular graph index exists in a specific branch
  bool _molIndexExistsInBranch(
    AtomIndex molIndex,
    TreeVertexIndex treeIndex
  ) const;

/* Major helper functions */
  /*! Adds missing adjacents to vertex, returning new and existing tree children
   *
   * Expands a tree node by checking the existing node children plus its parent
   * and comparing this against the adjacent indices from the molecule. Returns
   * all new and existing child tree vertex indices of the newly expanded node.
   * Pre-expansion existing children may come about due to the addition of
   * multiple bond order duplicate atoms.
   */
  std::vector<TreeVertexIndex> _expand(
    const TreeVertexIndex& index,
    const std::unordered_set<AtomIndex>& molIndicesInBranch
  );

  /*!
   * Since multiset's operator < does not actually USE the custom comparator
   * when comparing the contained values, we have to call
   * lexicographical_compare with a supplied comparator!
   */
  template<typename SetValueType, typename ComparatorType>
  bool _multisetCompare(
    const std::multiset<SetValueType, ComparatorType>& a,
    const std::multiset<SetValueType, ComparatorType>& b
  ) const {
    return std::lexicographical_compare(
      std::begin(a),
      std::end(a),
      std::begin(b),
      std::end(b),
      ComparatorType {*this}
    );
  }

  /*!
   * Compares all possible pairs of multisets of yet undecided branches and
   * adds less-than relationships to the orderingHelper in case new
   * relationships are discovered.
   *
   * NOTE: This is a common pattern in the application of sequence rules,
   * although the multiset value type and the applied comparator may vary, so
   * it is extracted here.
   *
   * @tparam SetValueType The multisets may contain different value types
   *   depending on the sequence rule being applied, but the semantics of how
   *   they are compared are always the same, so we abstract over this type.
   * @tparam ComparatorType The multisets may have differing comparators
   *   depending on the sequence rule being applied, but the semantics of how
   *   they are compared are always the same, so we abstract over this type.
   *
   * @param comparisonSets A mapping from branch indices that are ranked in the
   *   OrderDiscoveryHelper to multisets containing values that establish
   *   relative priority in the sequence rule currently being applied
   * @param undecidedSets This is a const reference to an up-to-date result
   *   value of the OrderDiscoveryHelper's getUndecidedSets() function. This
   *   is useful since undecidedSets is a commonly required variable in the
   *   application of sequence rules and can therefore be used across all of
   *   them
   * @param orderingHelper This is a reference to the OrderDiscoveryHelper
   *   currrently in use to rank branches of the acyclic tree. This function
   *   adds new relationships (if discovered) to this reference.
   */
  template<typename SetValueType, typename ComparatorType>
  void _compareBFSSets(
    const std::map<
      TreeVertexIndex,
      std::multiset<SetValueType, ComparatorType>
    >& comparisonSets,
    const std::vector<
      std::vector<TreeVertexIndex>
    >& undecidedSets,
    OrderDiscoveryHelper<TreeVertexIndex>& orderingHelper
  ) const {
    for(const auto& undecidedSet : undecidedSets) {
      temple::forEach(
        temple::adaptors::allPairs(undecidedSet),
        [&](const TreeVertexIndex a, const TreeVertexIndex b) {
          if(_multisetCompare(comparisonSets.at(a), comparisonSets.at(b))) {
            orderingHelper.addLessThanRelationship(a, b);
          } else if(_multisetCompare(comparisonSets.at(b), comparisonSets.at(a))) {
            orderingHelper.addLessThanRelationship(b, a);
          }
        }
      );
    }
  }

  std::string toString(TreeVertexIndex vertex) const;
  std::string toString(const TreeEdgeIndex& edge) const;
  std::string toString(const VariantType& vertexOrEdge) const;

  /* Some helper structs to get _runBFS to compile smoothly. Although the
   * boolean template parameters to _runBFS are logically constexpr, and any
   * false branches of if (boolean template parameter) are removed during
   * optimization, they must compile! This requirement is removed in C++17
   * if-constexpr structures.
   *
   * So EdgeInserter and VertexInserter's primary template (with boolean
   * parameter specifying whether they should perform a task or not) does
   * nothing, and there is a partial template specialization for the case
   * that the boolean template parameter is true, in which case the task is
   * performed.
   *
   * Classes are necessary since partial function specialization is impossible.
   *
   * When modernizing the code to C++17, remove those classes and employ much
   * more legible if-constexpr in all ifs in _runBFS using only constexpr
   * template parameters.
   */

  //! Helper struct to insert edges into the comparison set
  template<
    bool insertEdges,
    typename MultisetValueType,
    typename MultisetComparatorType
  > struct EdgeInserter {
    static void execute(
      std::map<
        TreeVertexIndex,
        std::multiset<MultisetValueType, MultisetComparatorType>
      >& /* comparisonSets */,
      const TreeVertexIndex& /* branchIndex */,
      const TreeEdgeIndex& /* edgeIndex */
    ) { /* do nothing */ }
  };

  //! Partial specialization that actually does the insertion
  template<
    typename MultisetValueType,
    typename MultisetComparatorType
  > struct EdgeInserter<true, MultisetValueType, MultisetComparatorType> {
    static void execute(
      std::map<
        TreeVertexIndex,
        std::multiset<MultisetValueType, MultisetComparatorType>
      >& comparisonSets,
      const TreeVertexIndex& branchIndex,
      const TreeEdgeIndex& edgeIndex
    ) {
      comparisonSets.at(branchIndex).insert(edgeIndex);
    }
  };

  //! Helper struct to insert vertices into the comparion set
  template<
    bool insertVertices,
    typename MultisetValueType,
    typename MultisetComparatorType
  > struct VertexInserter {
    static void execute(
      std::map<
        TreeVertexIndex,
        std::multiset<MultisetValueType, MultisetComparatorType>
      >& /* comparisonSets */,
      const TreeVertexIndex& /* branchIndex */,
      const TreeVertexIndex& /* vertexIndex */
    ) { /* do nothing */ }
  };

  //! Partial specialization that actually does the insertion
  template<
    typename MultisetValueType,
    typename MultisetComparatorType
  > struct VertexInserter<true, MultisetValueType, MultisetComparatorType> {
    static void execute(
      std::map<
        TreeVertexIndex,
        std::multiset<MultisetValueType, MultisetComparatorType>
      >& comparisonSets,
      const TreeVertexIndex& branchIndex,
      const TreeVertexIndex& vertexIndex
    ) {
      comparisonSets.at(branchIndex).insert(vertexIndex);
    }
  };

  /*! Performs a BFS traversal through the tree, performing ranking of a single rule
   *
   * Performs a BFS traversal through the tree, adding encountered edges and/or
   * vertices into a multiset with a custom comparator, and comparing those
   * multisets at every iteration to discover ordering relations between the
   * corresponding original branches.
   *
   * @tparam ruleNumber For logging purposes, indicate which sequence rule is
   *   being tested with this BFS traversal
   * @tparam BFSDownOnly If set true, only out-edges of any BFS seeds are
   *   considered for the multisets and seed continuations. This has the effect
   *   that the traversal is unidirectional, going "down" from the root of the
   *   tree, which is drawn at the top of the graphical representation
   *   (see dumpGraphviz). If set false, in- and out-edges are considered for
   *   BFS continuations, so that BFS proceeds in all directions simultaneously
   *   from every seed starting at the source index.
   * @tparam insertEdges If set true, edge indices are inserted into the
   *   multiset.
   * @tparam insertVertices If set true, vertex indices are inserted into the
   *   multiset.
   * @tparam MultisetValueType The value type of the multiset used for
   *   comparison.
   * @tparam MultisetComparatorType The Comparator supplied as constructing
   *   argument for the multiset and any lexicographical comparisons involving
   *   multiple multisets. Here, the edge or vertex properties being compared
   *   in the current sequence rule must be encoded.
   */
  template<
    unsigned ruleNumber,
    bool BFSDownOnly,
    bool insertEdges,
    bool insertVertices,
    typename MultisetValueType,
    typename MultisetComparatorType
  > void _runBFS(
    TreeVertexIndex sourceIndex,
    OrderDiscoveryHelper<TreeVertexIndex>& orderingHelper,
    const boost::optional<unsigned>& depthLimitOptional = boost::none
  ) const {
    /* Static safety checks */
    static_assert(
      insertEdges || insertVertices,
      "During BFS, you must insert either edges, vertices, or both into the "
      " comparison sets, but definitely not neither"
    );

    static_assert(
      (
        insertEdges
        && !std::is_same<MultisetValueType, TreeVertexIndex>::value
      ) || (
        insertVertices
        && !std::is_same<MultisetValueType, TreeEdgeIndex>::value
      ),
      "Multiset value type and insert booleans are mismatched"
    );

    // Check if depthLimit prohibits first comparison too
    if(depthLimitOptional.value_or(std::numeric_limits<unsigned>::max()) == 0) {
      return;
    }

    /* Declarations */

    /* For each branch to compare, keep a multiset of all the encountered
     * indices in the branch. The multiset keeps a descending order of the
     * indices according to sequence rule 1 and can be lexicographically
     * compared against one another using the SequenceRuleOneVertexComparator
     *
     * NOTE: Do not use operator [] to access members in comparisonSets. This
     * old form of accessible and assignable proxy object forces a
     * default-instantiation of the Comparator, which the
     * MultisetComparatorTypes cannot do. Always use .at().
     */
    std::map<
      TreeVertexIndex,
      std::multiset<MultisetValueType, MultisetComparatorType>
    > comparisonSets;

    /* For each branch to compare, keep a set of seed indices to follow in an
     * iteration. These are expanded in each iteration as long as the branch
     * they contain remains relevant (i.e. its relation to the other branches
     * is still unclear): They themselves are placed in the comparisonSet,
     * while their respective adjacents are the new seeds for the next
     * iteration.
     */
    std::map<
      TreeVertexIndex,
      std::vector<TreeVertexIndex>
    > seeds;

    /* In case the BFS in not down-only, we have to track which indices we
     * have visited to ensure BFS doesn't backtrack or reuse vertices.
     */
    std::unordered_set<TreeVertexIndex> visitedVertices;


    /* Initialization */
    EdgeInserter<insertEdges, MultisetValueType, MultisetComparatorType> edgeInserter;
    VertexInserter<insertVertices, MultisetValueType, MultisetComparatorType> vertexInserter;

    if(!BFSDownOnly) { // Initialization is only necessary in this case
      visitedVertices.insert(sourceIndex);
    }

    auto undecidedSets = orderingHelper.getUndecidedSets();

    for(const auto& undecidedSet : undecidedSets) {
      for(const auto& undecidedBranch : undecidedSet) {
        comparisonSets.emplace(
          undecidedBranch,
          *this
        );

        vertexInserter.execute(
          comparisonSets,
          undecidedBranch,
          undecidedBranch
        );

        if(insertEdges) {
          if(BFSDownOnly) {
            edgeInserter.execute(
              comparisonSets,
              undecidedBranch,
              boost::edge(sourceIndex, undecidedBranch, _tree).first
            );
          } else {
            // Need to be direction-agnostic
            auto forwardEdge = boost::edge(sourceIndex, undecidedBranch, _tree);
            auto backwardEdge = boost::edge(undecidedBranch, sourceIndex, _tree);

            edgeInserter.execute(
              comparisonSets,
              undecidedBranch,
              (
                forwardEdge.second
                ? forwardEdge.first
                : backwardEdge.first
              )
            );
          }
        }

        seeds[undecidedBranch] = {undecidedBranch};
      }
    }

    // Compare undecided multisets, and add any discoveries to the ordering helper
    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    // Update the undecided sets
    undecidedSets = orderingHelper.getUndecidedSets();

    if /* C++17 constexpr */ (buildTypeIsDebug) {
      // Write debug graph files if the corresponding log particular is set
      if(
        Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0
        && _notEmpty(comparisonSets)
      ) {
        std::string header = (
          (
            BFSDownOnly
            ? ""s
            : "aux "s
          ) + "R"s + std::to_string(ruleNumber)
        );

        _writeGraphvizFiles({
          _adaptedMolGraphviz,
          dumpGraphviz(
            header,
            {sourceIndex},
            _collectSeeds(seeds, undecidedSets)
          ),
          _makeBFSStateGraph(
            header,
            sourceIndex,
            comparisonSets,
            undecidedSets
          ),
          orderingHelper.dumpGraphviz()
        });
      }
    }

    unsigned depth = 1;

    /* Loop BFS */

    /* As long as there are undecided sets and seeds whose expansion could be
     * relevant for those undecided sets, continue BFS
     */
    while(
      // Undecided indices remain
      !undecidedSets.empty()
      // Seeds exist that are relevant to the undecided sets
      && _relevantSeeds(seeds, undecidedSets)
      && depth < depthLimitOptional.value_or(std::numeric_limits<unsigned>::max())
    ) {
      // BFS Step
      for(const auto& undecidedSet: undecidedSets) {
        for(const auto& undecidedBranch: undecidedSet) {
          /* Clear the comparison set at undecided branches prior to inserting
           * anything. This doesn't change anything semantically, since
           * everything from the prior iteration (or initialization) must have
           * compared equal within this undecidedSet.
           *
           * Makes for cleaner visualization and reduces the number of
           * Comparator calls during insert and lexicographical_compare.
           */
          comparisonSets.at(undecidedBranch).clear();

          std::vector<TreeVertexIndex> newSeeds;

          for(const auto& seed : seeds.at(undecidedBranch)) {

            if(!BFSDownOnly) {
              // Mark as visited
              visitedVertices.insert(seed);

              // In-edges are only relevant for omnidirectional BFS
              for(
                auto inIterPair = boost::in_edges(seed, _tree);
                inIterPair.first != inIterPair.second;
                ++inIterPair.first
              ) {
                const auto& inEdge = *inIterPair.first;

                auto edgeSource = boost::source(inEdge, _tree);

                // Check if already placed
                if(visitedVertices.count(edgeSource) == 0) {
                  vertexInserter.execute(
                    comparisonSets,
                    undecidedBranch,
                    edgeSource
                  );

                  edgeInserter.execute(
                    comparisonSets,
                    undecidedBranch,
                    inEdge
                  );

                  newSeeds.push_back(edgeSource);
                }
              }
            }

            // Out edges are considered in all cases
            for( // Out-edges
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdge = *outIterPair.first;

              auto edgeTarget = boost::target(outEdge, _tree);

              // Skip this vertex if in omnidirectional BFS and already-seen node
              if(!BFSDownOnly && visitedVertices.count(edgeTarget) > 0) {
                continue;
              }

              vertexInserter.execute(
                comparisonSets,
                undecidedBranch,
                edgeTarget
              );

              edgeInserter.execute(
                comparisonSets,
                undecidedBranch,
                outEdge
              );

              // Add out edge target to seeds only if non-terminal
              if(!_tree[edgeTarget].isDuplicate) {
                newSeeds.push_back(edgeTarget);
              }
            }
          }

          // Overwrite seeds
          seeds.at(undecidedBranch) = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

      // Increment depth
      ++depth;

      if /* C++17 constexpr */ (buildTypeIsDebug) {
        if(
          Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0
          && _notEmpty(comparisonSets)
        ) {
          std::string header;
          if(!BFSDownOnly) {
            header += "aux "s;
          }
          header += "R"s + std::to_string(ruleNumber);

          _writeGraphvizFiles({
            _adaptedMolGraphviz,
            dumpGraphviz(
              header,
              {sourceIndex},
              _collectSeeds(seeds, undecidedSets)
            ),
            _makeBFSStateGraph(
              header,
              sourceIndex,
              comparisonSets,
              undecidedSets
            ),
            orderingHelper.dumpGraphviz()
          });
        }
      }
    }
  }

  // Flatten the undecided branches' seeds into a single set
  static std::unordered_set<TreeVertexIndex> _collectSeeds(
    const std::map<
      TreeVertexIndex,
      std::vector<TreeVertexIndex>
    >& seeds,
    const std::vector<
      std::vector<TreeVertexIndex>
    >& undecidedSets
  );

  //! Creates a graphviz representation of the BFS state
  template<typename ValueType, typename ComparatorType>
  std::string _makeBFSStateGraph(
    const std::string& title,
    const TreeVertexIndex& base,
    const std::map<
      TreeVertexIndex,
      std::multiset<ValueType, ComparatorType>
    >& comparisonSet,
    const std::vector<
      std::vector<TreeVertexIndex>
    >& undecidedSets
  ) const {

    std::set<TreeVertexIndex> activeBranches;
    for(const auto& set : undecidedSets) {
      for(const auto& value : set) {
        activeBranches.insert(value);
      }
    }

    struct DisplayGraphVertexData {
      std::string representation;
    };

    using DisplayGraph = boost::adjacency_list<
      /* OutEdgeListS = Type of Container for edges of a vertex
       * Options: vector, list, slist, set, multiset, unordered_set
       * Choice: setS
       */
      boost::setS,
      /* VertexListS = Type of Container for vertices
       * Options: vector, list, slist, set, multiset, unordered_set
       * Choice: vecS, removing vertices does not occur
       */
      boost::vecS,
      /* DirectedS = Is the graph directed or not?
       * Choice: Directed, do not need bidirectional
       */
      boost::directedS,
      /* VertexProperty = What information is stored about vertices?
       * Choice: Atom, containing an index and an element type
       */
      DisplayGraphVertexData
    >;

    std::set<TreeVertexIndex> colorVertices;
    std::set<TreeVertexIndex> squareVertices;

    DisplayGraph graph;

    auto baseVertex = boost::add_vertex(graph);
    graph[baseVertex].representation = toString(base);

    squareVertices.insert(baseVertex);

    for(const auto& branchIterPair : comparisonSet) {
      const auto& treeBranchIndex = branchIterPair.first;
      const auto& treeBranchMultiset = branchIterPair.second;

      auto branchVertex = boost::add_vertex(graph);
      graph[branchVertex].representation = toString(treeBranchIndex);

      squareVertices.insert(branchVertex);

      if(activeBranches.count(treeBranchIndex) > 0) {
        colorVertices.insert(branchVertex);
      }

      boost::add_edge(baseVertex, branchVertex, graph);

      typename DisplayGraph::vertex_descriptor continuation = branchVertex;
      for(const auto& multisetValue : treeBranchMultiset) {
        auto newVertex = boost::add_vertex(graph);
        graph[newVertex].representation = toString(multisetValue);

        boost::add_edge(continuation, newVertex, graph);
        continuation = newVertex;
      }
    }

    class SmallGraphWriter {
    private:
      const DisplayGraph& _graphBase;
      const std::string _title;
      const std::set<TreeVertexIndex> _colorVertices;
      const std::set<TreeVertexIndex> _squareVertices;

    public:
      SmallGraphWriter(
        const DisplayGraph& base,
        std::string title,
        std::set<TreeVertexIndex> colorVertices,
        std::set<TreeVertexIndex> squareVertices
      ) : _graphBase(base),
          _title(std::move(title)),
          _colorVertices(std::move(colorVertices)),
          _squareVertices(std::move(squareVertices))
      {}

      void operator() (std::ostream& os) const {
        os << "  graph [fontname = \"Arial\", layout=\"dot\"];\n"
          << "  node [fontname = \"Arial\", shape = circle, style = filled];\n"
          << "  edge [fontname = \"Arial\"];\n"
          << R"(  labelloc="t"; label=")" << _title << "\"\n";
      }

      void operator() (
        std::ostream& os,
        const typename DisplayGraph::vertex_descriptor& vertexIndex
      ) const {
        if(vertexIndex == rootIndex) {
          os << R"([label=")" << _title << "\\n-\\n" << _graphBase[rootIndex].representation << R"(")";
        } else {
          os << R"([label=")" << _graphBase[vertexIndex].representation << R"(")";
        }

        if(_colorVertices.count(vertexIndex) > 0) {
          os << R"(, fillcolor="tomato", fontcolor="white")";
        }

        if(_squareVertices.count(vertexIndex) > 0) {
          os << R"(, shape="square")";
        }

        os << R"(])";
      }

      void operator() (
        std::ostream& /* os */,
        const typename DisplayGraph::edge_descriptor& /* edge */
      ) const { /* do nothing particular for edges */ }

    };

    SmallGraphWriter propertyWriter {graph, title, colorVertices, squareVertices};

    std::stringstream ss;

    boost::write_graphviz(
      ss,
      graph,
      propertyWriter,
      propertyWriter,
      propertyWriter
    );

    return ss.str();
  }

  //! Creates graphviz representation of like/unlike pairing in sequence rule 4B
  std::string _make4BGraph(
    const TreeVertexIndex& sourceIndex,
    const std::map<
      TreeVertexIndex,
      std::set<VariantType>
    >& representativeStereodescriptors,
    const TreeVertexIndex& branchA,
    const TreeVertexIndex& branchB,
    const std::vector<
      std::vector<VariantType>
    >& branchAOrders,
    const std::vector<
      std::vector<VariantType>
    >& branchBOrders,
    const std::vector<
      std::vector<VariantType>
    >::reverse_iterator& branchAIter,
    const std::vector<
      std::vector<VariantType>
    >::reverse_iterator& branchBIter
  ) const;

  //! Maps sets returned by an OrderDiscoveryHelper from tree indices to atom indices
  std::vector<
    std::vector<AtomIndex>
  > _mapToAtomIndices(
    const std::vector<
      std::vector<TreeVertexIndex>
    >& treeRankingSets
  ) const;

  /*!
   * In all BFS-like iterations, we need to check that there are suitable
   * seeds to continue the BFS expansion for all branches that are yet
   * undifferentiated. This is required in pretty much every sequence rule
   * application, so it is generalized here.
   *
   * @param seeds A mapping from branch indices that are being ranked to
   *   a set of tree vertices that are suitable continuations of the current
   *   BFS iteration.
   * @param undecidedSets This is a const reference to an up-to-date result
   *   value of the OrderDiscoveryHelper's getUndecidedSets() function. This
   *   is useful since undecidedSets is a commonly required variable in the
   *   application of sequence rules and can therefore be used across all of
   *   them
   * @returns Whether any of the branches that are yet considered equal (and
   *   thus together in one of the undecided sets) have any seeds for another
   *   BFS iteration
   */
  static bool _relevantSeeds(
    const std::map<
      TreeVertexIndex,
      std::vector<TreeVertexIndex>
    >& seeds,
    const std::vector<
      std::vector<TreeVertexIndex>
    >& undecidedSets
  );

  /*! Ranks the direct substituents of the root atom, applying sequence rules 2-5
   *
   * This function ranks the direct substituents of the atom this RankingTree
   * was instantiated upon by the sequential application of the 2013 IUPAC Blue
   * Book sequence rules 2-5 (the ctor applies #1). They are somewhat adapted
   * since the priority of asymmetric centers of higher (and also lower)
   * shape must also be considered because transition metal chemistry is
   * also included in this library.
   */
  void _applySequenceRules(
    const boost::optional<AngstromWrapper>& positionsOption
  );

  /*!
   * This function ranks the direct substituents of a selected central tree
   * vertex by the sequential application of the 2013 IUPAC Blue Book sequence
   * rules. They are somewhat adapted since the priority of asymmetric centers
   * of higher (and also lower) shape must also be considered because
   * transition metal chemistry is also included in this library.
   *
   * It returns a sorted vector of vectors, in which every sub-vector
   * represents a set of equal-priority tree vertices. The sorting is ascending,
   * meaning from lowest priority to highest priority.
   *
   * The main differences to rank() are that all sequence rule BFS progressions
   * can begin at any vertex in the tree and therefore must consider both in-
   * and out-edges at every vertex. Additionally, this function assumes that
   * any auxiliary stereodescriptors have already been instantiated, while
   * rank() explicitly instantiates them prior to the application of sequence
   * rule three.
   *
   * @param sourceIndex The tree vertex from which to rank substituents
   * @param adjacentsToRank A list of adjacent tree vertices to rank
   * @param depthLimitOptional An optional limitation on depth of sequence rule
   *   application
   *
   * @returns A ranked nested list of tree vertex sets as a ragged vector
   */
  std::vector<
    std::vector<TreeVertexIndex>
  > _auxiliaryApplySequenceRules(
    TreeVertexIndex sourceIndex,
    const std::vector<TreeVertexIndex>& adjacentsToRank,
    const boost::optional<unsigned>& depthLimitOptional = boost::none
  ) const;

public:
//!@name Special member functions
//!@{
  /*! @brief Performs ranking of a central atom's substituents
   *
   * @complexity{Theoretical complexity is unclear to me. No idea if the IUPAC
   * sequence rules imply an upper bound on theoretical complexity. Given a
   * sufficiently complicated molecule, a single ranking can be immensely
   * complicated and time intensive. For practical cases however, ranking is
   * constant-time.}
   */
  RankingTree(
    const OuterGraph& graph,
    const Cycles& cycles,
    const StereopermutatorList& stereopermutators,
    std::string molGraphviz,
    AtomIndex atomToRank,
    const std::vector<AtomIndex>& excludeIndices = {},
    ExpansionOption expansionMethod = ExpansionOption::OnlyRequiredBranches,
    const boost::optional<AngstromWrapper>& positionsOption = boost::none
  );
//!@}

//!@name Information
//!@{
  /*! Fetches the ranked result
   *
   * Returns the ranked result as a sorted vector of vectors, in which every
   * sub-vector represents a set of equal-priority substituents. The sorting is
   * ascending, meaning from lowest priority to highest priority.
   */
  std::vector<
    std::vector<AtomIndex>
  > getRanked() const;

  /*! Returns an annotated graphviz graph of the tree
   *
   * Creates a graphviz representation of the tree, with optional title string,
   * and sets of vertices represented as squares or colored in.
   */
  std::string dumpGraphviz(
    const std::string& title = "",
    const std::unordered_set<TreeVertexIndex>& squareVertices = {},
    const std::unordered_set<TreeVertexIndex>& colorVertices = {},
    const std::set<TreeEdgeIndex>& colorEdges = {}
  ) const;
//!@}
};

} // namespace molassembler

} // namespace Scine

#endif
