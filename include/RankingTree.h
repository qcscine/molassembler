#ifndef INCLUDE_MOLECULE_MANIP_RANKING_HIERARCHICAL_TREE_H
#define INCLUDE_MOLECULE_MANIP_RANKING_HIERARCHICAL_TREE_H

#include "BondDistance.h"
#include "Log.h"
#include "Molecule.h"
#include "OrderDiscoveryHelper.h"

/* TODO
 * - Pseudo-asymmetry considerations
 *   - Is the propagation of pseudoasymmetry REALLY necessary? It makes sense
 *     to permute a stereocenter designated as pseudo-asymmetric. Dunno if
 *     creating the exact diastereomer of a molecule will ever be needed as an
 *     operation (if so, then yes, pseudoasymmetry pop is needed). Creating
 *     ALL stereomers of a molecule can be done without pseudoasymmetry.
 *     All other points below are only relevant in case pseudo-asymmetry needs
 *     to be included.
 *   - Stereocenter interface change to support pseudo-asymmetry tag (?)
 *   - Ranking function interface change to propagate pseudo-asymmetry result
 * - Optimizations / Refactors
 *   - 4A could be refactored into _runBFS, if 4B could collect the variants it
 *     needs to rank itself. Don't know if an improvement to avoid BFS-ing again
 *     to find stereocenters over potential gains from clearing comparisonSets
 *     every BFS step
 * - Instantiation of EZStereocenters on edge does not keep CNStereocenters
 *   from being instantiated above the edge in sequence rule 3 prep
 *   (see 2Z... file ranking), probably innocuous, but unnecessary
 * - OrderDiscoveryHelper function naming may be inconsistent. Does
 *   addLessThanRelationship really semantically do what it says?
 *   Loads of inverted comparisons here where it doesn't make sense to me
 *   anymore, in particular at 4B in counting like pairs at every position
 * - TODOs strewn about
 */

/*! @file
 *
 * Centerpoint of library ranking algorithm. Implements the RankingTree class,
 * which can be instantiated on any atomic index in a Molecule, which then
 * splits the Molecule into an acyclic tree, which can then be used to rank
 * that atom's direct substituents according to IUPAC-like sequence rules.
 */

namespace MoleculeManip {

/*!
 * Central class for unified IUPAC-like ranking of organic and inorganic
 * structures. The general procedure is that any cycles are broken down in a BFS
 * iteration through the molecular graph at the atom vertex of instantiation. 
 * Afterwards, all branches are expanded in a DFS-like manner to ensure that
 * cycle closures are properly reported as duplicate atoms. In the rank() member
 * function, the actual ranking as demanded by the IUPAC sequence rules (adapted
 * for the mixed inorganic/organic molecular graphs) is performed.
 */
class RankingTree {
#ifndef NEDBUG
private:
  static unsigned _debugMessageCounter;

  static std::string _adaptMolGraph(std::string molGraph);

  static void _writeGraphvizFiles(
    const std::vector<std::string>& graphvizStrings
  );
#endif

public:
  enum ExpansionOption {
    Optimized,
    Full
  };

  //! Data class that sets which supplementary data is stored for a tree vertex
  struct VertexData {
    AtomIndexType molIndex;
    bool isDuplicate;

    boost::optional<Stereocenters::CNStereocenter> stereocenterOption;

    // TODO not quite ready for this yet
    // boost::optional<double> deviantAtomicNumber;
  };

  //! Data class that sets which supplementary data is stored for a tree edge
  struct EdgeData {
    boost::optional<Stereocenters::EZStereocenter> stereocenterOption;
  };

  //! The BGL Graph type used to store the tree
  using TreeGraphType = boost::adjacency_list<
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

private:
  //! Type of tree vertex accessor
  using TreeVertexIndex = TreeGraphType::vertex_descriptor;
  //! Type of tree edge accessor
  using TreeEdgeIndex = TreeGraphType::edge_descriptor;
  //! Variant type of both
  using VariantType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

  static constexpr TreeVertexIndex rootIndex = 0;

/* State */
  //! The BGL Graph representing the acyclic tree
  TreeGraphType _tree;
  /*!
   * A set of edge indices gathered construction for all tree edges representing
   * a molecular bond of order two.
   */
  std::set<TreeEdgeIndex> _doubleBondEdges;
  //! The helper instance for discovering the ordering of the to-rank branches
  OrderDiscoveryHelper<TreeVertexIndex> _branchOrderingHelper;

  //! Overall order discoveries from _auxiliary calls
  OrderDiscoveryHelper<TreeVertexIndex> _allOrdering;

  // Closures
  const Molecule& _moleculeRef;

/* Minor helper classes and functions */
  //! Helper class to write a graphviz representation of the generated tree
  class GraphvizWriter;

  //! Returns the parent of a node. Fails if called on the root!
  TreeVertexIndex _parent(const TreeVertexIndex& index) const;

  //! Returns the direct descendants of a tree node
  std::set<TreeVertexIndex> _children(const TreeVertexIndex& index) const;

  //! Returns a set of all adjacent vertices (in- and out-adjacents)
  std::set<TreeVertexIndex> _adjacents(const TreeVertexIndex& index) const;

  //! Returns a set of all adjacent edges (in- and out-edges)
  std::set<TreeEdgeIndex> _adjacentEdges(const TreeVertexIndex& index) const;

  unsigned _adjacentTerminalHydrogens(const TreeVertexIndex& index) const;

  bool _isBondSplitDuplicateVertex(const TreeVertexIndex& index) const;

  bool _isCycleClosureDuplicateVertex(const TreeVertexIndex& index) const;

  std::set<TreeVertexIndex> _auxiliaryAdjacentsToRank(
    const TreeVertexIndex& sourceIndex,
    const std::set<TreeVertexIndex>& excludedIndices
  ) const;

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
  std::set<TreeVertexIndex> _treeIndicesInBranch(TreeVertexIndex index) const;

  //! Returns all mol indices in the branch from the specified index up to root
  std::set<AtomIndexType> _molIndicesInBranch(const TreeVertexIndex& index) const;

  //! Returns a depth measure of a duplicate vertex for sequence rule 1
  unsigned _duplicateDepth(TreeVertexIndex index) const;

  //! Returns the depth of a node in the tree
  unsigned _depthOfNode(TreeVertexIndex index) const;

  //! Returns a mixed depth measure for ranking both vertices and edges
  unsigned _mixedDepth(const TreeVertexIndex& vertexIndex) const;

  //! Returns a mixed depth measure for ranking both vertices and edges
  unsigned _mixedDepth(const TreeEdgeIndex& edgeIndex) const;

  //! A data class, stores both the junction vertex and paths from the source vertices.
  struct JunctionInfo;

  /*!
   * Returns the deepest vertex in the tree in whose child branches both a and
   * b are located
   */
  JunctionInfo _junction(const TreeVertexIndex& a, const TreeVertexIndex& b) const;

  //! Returns whether a molecular graph index exists in a specific branch
  bool _molIndexExistsInBranch(
    const AtomIndexType& molIndex,
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
    const std::set<AtomIndexType>& molIndicesInBranch
  );

  /*!
   * Since multiset's operator < does not actually USE the custom comparator
   * when comparing the contained values, we have to call
   * lexicographical_compare, supplying the comparator!
   */
  template<typename SetValueType, typename ComparatorType>
  bool _multisetCompare(
    const std::multiset<SetValueType, ComparatorType>& a,
    const std::multiset<SetValueType, ComparatorType>& b
  ) const {
    return std::lexicographical_compare(
      a.begin(),
      a.end(),
      b.begin(),
      b.end(),
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
      TemplateMagic::forAllPairs(
        undecidedSet,
        [&](const TreeVertexIndex& a, const TreeVertexIndex& b) {
          if(_multisetCompare(comparisonSets.at(a), comparisonSets.at(b))) {
            orderingHelper.addLessThanRelationship(a, b);
          } else if(_multisetCompare(comparisonSets.at(b), comparisonSets.at(a))) {
            orderingHelper.addLessThanRelationship(b, a);
          }
        }
      );
    }
  }

  template<typename ValueType>
  std::string toString(const ValueType& value) const {
    return ""s;
  }

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
      >& comparisonSets,
      const TreeVertexIndex& branchIndex,
      const TreeEdgeIndex& edgeIndex
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
      >& comparisonSets,
      const TreeVertexIndex& branchIndex,
      const TreeVertexIndex& vertexIndex
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

  /* When modernizing to C++17, see the comments regarding EdgeInserter and
   * VertexInserter above
   */
  /*! Performs a BFS traversal through the tree
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
    OrderDiscoveryHelper<TreeVertexIndex>& orderingHelper
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
    std::set<TreeVertexIndex> visitedVertices;


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
        /*if(insertVertices) {
          comparisonSets.at(undecidedBranch).insert(undecidedBranch);
        }*/

        if(insertEdges) {
          if(BFSDownOnly) {
            edgeInserter.execute(
              comparisonSets,
              undecidedBranch,
              boost::edge(sourceIndex, undecidedBranch, _tree).first
            );

            /*comparisonSets.at(undecidedBranch).insert(
              boost::edge(sourceIndex, undecidedBranch, _tree).first
            );*/
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

            /*comparisonSets.at(undecidedBranch).insert(
              forwardEdge.second
              ? forwardEdge.first
              : backwardEdge.first
            );*/
          }
        }

        seeds[undecidedBranch] = {undecidedBranch};
      }
    }

    // Compare undecided multisets, and add any discoveries to the ordering helper
    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    // Update the undecided sets
    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    // Write debug graph files if the corresponding log particular is set
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      std::string header = (
        (
          BFSDownOnly
          ? ""s
          : "aux "s
        ) + "R"s + std::to_string(ruleNumber)
      );

      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz(
          header,
          {sourceIndex},
          _collectSeeds(seeds, undecidedSets)
        ),
        _makeGraph(
          header,
          sourceIndex,
          comparisonSets,
          undecidedSets
        ),
        orderingHelper.dumpGraphviz()
      });
    }
#endif


    /* Loop BFS */

    /* As long as there are undecided sets and seeds whose expansion could be
     * relevant for those undecided sets, continue BFS
     */
    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
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
                  /*if(insertVertices) {
                    comparisonSets.at(undecidedBranch).emplace(edgeSource);
                  }*/

                  edgeInserter.execute(
                    comparisonSets,
                    undecidedBranch,
                    inEdge
                  );
                  /*if(insertEdges) {
                    comparisonSets.at(undecidedBranch).emplace(inEdge);
                  }*/

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
              /*if(insertVertices) {
                comparisonSets.at(undecidedBranch).emplace(edgeTarget);
              }*/

              edgeInserter.execute(
                comparisonSets,
                undecidedBranch,
                outEdge
              );
              /*if(insertEdges) {
                comparisonSets.at(undecidedBranch).emplace(outEdge);
              }*/

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

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        std::string header;
        if(!BFSDownOnly) {
          header += "aux "s;
        }
        header += "R"s + std::to_string(ruleNumber);

        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz(
            header,
            {sourceIndex},
            _collectSeeds(seeds, undecidedSets)
          ),
          _makeGraph(
            header,
            sourceIndex,
            comparisonSets,
            undecidedSets
          ),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }
  }

  // Flatten the undecided branches' seeds into a single set
  static std::set<TreeVertexIndex> _collectSeeds(
    const std::map<
      TreeVertexIndex,
      std::vector<TreeVertexIndex>
    >& seeds,
    const std::vector<
      std::vector<TreeVertexIndex>
    >& undecidedSets
  );

  template<typename ValueType, typename ComparatorType>
  std::string _makeGraph(
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
        const std::string& title,
        const std::set<TreeVertexIndex> colorVertices,
        const std::set<TreeVertexIndex> squareVertices
      ) : _graphBase(base),
          _title(title),
          _colorVertices(colorVertices),
          _squareVertices(squareVertices)
      {}

      void operator() (std::ostream& os) const {
        os << "  graph [fontname = \"Arial\"];\n"
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
        std::ostream&,
        const typename DisplayGraph::edge_descriptor&
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

  std::vector<
    std::vector<AtomIndexType>
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
   * symmetry must also be considered because transition metal chemistry is
   * also included in this library.
   */
  void _applySequenceRules(
    const boost::optional<Delib::PositionCollection>& positionsOption
  );

  /*!
   * This function ranks the direct substituents of a selected central tree
   * vertex by the sequential application of the 2013 IUPAC Blue Book sequence
   * rules. They are somewhat adapted since the priority of asymmetric centers
   * of higher (and also lower) symmetry must also be considered because
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
   */
  std::vector<
    std::vector<TreeVertexIndex>
  > _auxiliaryApplySequenceRules(
    const TreeVertexIndex& sourceIndex,
    const std::set<TreeVertexIndex>& adjacentsToRank
  ) const;


public:
/* Sequence rule comparator classes for use in std::multiset */

  /*! 
   * Comparator class for comparing individual tree vertex indices according to
   * IUPAC Sequence rule one. Handles all combinations of duplicate and
   * non-duplicate nodes.
   *
   * NOTE: This is a non-default-instantiable Comparator, meaning care must be
   * taken in the instantiation of the STL Container using this to avoid move
   * and copy assignment operators. You must use in-place-construction!
   */
  struct SequenceRuleOneVertexComparator;

  /*! 
   * Comparator class for comparing individual tree edges according to IUPAC
   * sequence rule three.
   *
   * NOTE: This is a non-default-instantiable Comparator, meaning care must be
   * taken in the instantiation of the STL Container using this to avoid move
   * and copy assignment operators. You must use in-place-construction!
   */
  struct SequenceRuleThreeEdgeComparator;

  /*! 
   * Comparator class for comparing both tree vertices and edges according to
   * IUPAC sequence rule four.
   *
   * NOTE: This is a non-default-instantiable Comparator, meaning care must be
   * taken in the instantiation of the STL Container using this to avoid move
   * and copy assignment operators. You must use in-place-construction!
   */
  struct SequenceRuleFourVariantComparator;

  /*! 
   * Comparator class for comparing both tree vertices and edges according to
   * IUPAC sequence rule five.
   *
   * NOTE: This is a non-default-instantiable Comparator, meaning care must be
   * taken in the instantiation of the STL Container using this to avoid move
   * and copy assignment operators. You must use in-place-construction!
   */
  struct SequenceRuleFiveVariantComparator;


/* Variant visitor helper classes */

  //! Returns whether the vertex or edge has an instantiated Stereocenter
  class VariantHasInstantiatedStereocenter;

  //! Returns the mixed depth (see _mixedDepth functions) of the vertex or edge
  class VariantDepth;

  //! Returns the source node of a vertex (identity) or an edge (source node)
  class VariantSourceNode;

  //! Returns a string representation of the stereocenter on vertex or edge
  class VariantStereocenterStringRepresentation;

  /*!
   * Returns whether two variants' stereocenters form a like pair or not
   * according to the sub-rules in IUPAC sequence rule four.
   */
  class VariantLikePair;


/* Constructors */

  /*! Creates an acyclic tree from the molecular graph.
   *
   * Splits any cycles of the molecular graph into an acyclic tree with its root
   * at the atom instantiated upon.
   */
  RankingTree(
    const Molecule& molecule,
    const AtomIndexType& atomToRank,
    const std::set<AtomIndexType>& excludeIndices = {},
    const ExpansionOption& expansionMethod = ExpansionOption::Optimized,
    const boost::optional<Delib::PositionCollection>& positionsOption = boost::none
  );

  /*! Fetches the ranked result
   *
   * Returns the ranked result as a sorted vector of vectors, in which every
   * sub-vector represents a set of equal-priority substituents. The sorting is
   * ascending, meaning from lowest priority to highest priority.
   */
  std::vector<
    std::vector<AtomIndexType>
  > getRanked() const;

  /*! Returns an annotated graphviz graph of the tree
   *
   * Creates a graphviz representation of the tree, with optional title string,
   * and sets of vertices represented as squares or colored in.
   */
  std::string dumpGraphviz(
    const std::string& title = "",
    const std::set<TreeVertexIndex>& squareVertices = {},
    const std::set<TreeVertexIndex>& colorVertices = {},
    const std::set<TreeEdgeIndex>& colorEdges = {}
  ) const;

  const TreeGraphType& getGraph() const;
};

} // namespace MoleculeManip

#endif
