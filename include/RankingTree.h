#ifndef INCLUDE_MOLECULE_MANIP_RANKING_HIERARCHICAL_TREE_H
#define INCLUDE_MOLECULE_MANIP_RANKING_HIERARCHICAL_TREE_H

#include "BondDistance.h"
#include "Molecule.h"
#include "OrderDiscoveryHelper.h"

/* TODO
 * - Stereocenter interface change needed to consume ascending list of equal
 *   sets, plus modify this everywhere in the codebase, plus test correctness
 * - Stereocenter interface change to support pseudo-asymmetry tag
 * - Ranking function interface change to propagate pseudo-asymmetry result
 * - Seeds shouldn't be sets, but vectors (no need for set properties)
 * - Visualize sequence rule 4B
 * - Interface change to support using coordinates or prior instantiated
 *   stereocenters to assign auxiliary stereocenters
 * - Since most rankings are done after a few sequence rule 1 iterations, it may
 *   not be necessary to even build the entire tree from the molecule. This will
 *   be *absolutely necessary* for large molecules such as proteins or molecules
 *   with many fused cycles.
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

/* State */
  //! The BGL Graph representing the acyclic tree
  TreeGraphType _tree;
  /*!
   * A set of edge indices gathered construction for all tree edges representing
   * a molecular bond of order two.
   */
  std::set<TreeEdgeIndex> _doubleBondEdges;

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

  /*! 
   * Returns the node degree, excluding edges to or from nodes marked as
   * duplicate
   */
  unsigned _nonDuplicateDegree(const TreeVertexIndex& index) const;

  //! Adds the required duplicate atoms for bonds of high bond orders
  void _addBondOrderDuplicates(
    const TreeVertexIndex& treeSource,
    const TreeVertexIndex& treeTarget
  );

  //! Returns all indices in the branch from the specified index up to root
  std::set<AtomIndexType> _molIndicesInBranch(TreeVertexIndex index) const;

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
  /*! Expand a node onto its missing adjacencies by depth-first recursive
   * progression
   */
  void _DFSExpand(
    const TreeVertexIndex& index,
    const std::set<AtomIndexType>& molIndicesInBranch
  );

  /*!
   * Since multiset's operator < does not actually USE the custom comparator
   * when comparing the contained values, we have to call
   * lexicographical_compare ourselves, supplying the comparator!
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
        os << R"([label=")" << _graphBase[vertexIndex].representation << R"(")";

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
      std::set<TreeVertexIndex>
    >& seeds,
    const std::vector<
      std::vector<TreeVertexIndex>
    >& undecidedSets
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
  > _omniDirectionalRank(
    const TreeVertexIndex& sourceIndex,
    const std::set<TreeVertexIndex>& excludeIndices
  ) const;

public:
/* BFS Visitor for initial cycle splitting */
  /*!
   * An adaptation of GraphAlgorithms' TreeGenerator, which created a pointer-
   * based Tree (see Tree.h). This generates a tree within a BGL Bidirectional
   * graph, but the operating principle is identical.
   */
  class RankingTreeGenerator;


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
    const AtomIndexType& atomToRank
  );

  /*! Ranks the direct substituents of the root atom
   *
   * This function ranks the direct substituents of the atom this RankingTree
   * was instantiated upon by the sequential application of the 2013 IUPAC Blue
   * Book sequence rules. They are somewhat adapted since the priority of
   * asymmetric centers of higher (and also lower) symmetry must also be
   * considered because transition metal chemistry is also included in this
   * library.
   *
   * It returns a sorted vector of vectors, in which every sub-vector
   * represents a set of equal-priority substituents. The sorting is ascending,
   * meaning from lowest priority to highest priority.
   */
  std::vector<
    std::vector<AtomIndexType>
  > rank(const std::set<AtomIndexType>& excludeIndices = {});

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
