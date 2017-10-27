#include "RankingTree.h"

#include "boost/graph/breadth_first_search.hpp"
#include "boost/algorithm/string/replace.hpp"
#include "template_magic/Boost.h"

#include "Delib/ElementInfo.h"

#include "Log.h"
#include "MolGraphWriter.h"

#include <fstream>
#include <iostream>

namespace MoleculeManip {

/* Graph splitting class */
class RankingTree::RankingTreeGenerator : public boost::default_bfs_visitor {
public:
  using KeyToNodeListMap = std::map<
    AtomIndexType,
    std::vector<TreeVertexIndex>
  >;

private:
/* Internal state */
  //! Base reference to tree being generated
  RankingTree& _treeBaseRef;
  //! Map to keep track of depth from starting vertex
  std::map<AtomIndexType, unsigned> _depthMap;
  //! Pointer to Tree root on heap (where we generate it)
  TreeVertexIndex _rootIndex;
  //! Mapping of atomic index to nodes in the Tree currently being generated
  KeyToNodeListMap _existingNodePtrMap;
  //! List of edges to add in current depth
  std::vector<
    std::pair<AtomIndexType, AtomIndexType>
  > _edgesToAddBeforeDescent;


public:
  struct EarlyExit {};

  RankingTreeGenerator(
    RankingTree& baseTree,
    const AtomIndexType& startingFrom,
    const TreeVertexIndex& rootIndex
  ) : _treeBaseRef(baseTree),
      _depthMap({
        {startingFrom, 0}
      }),
      _rootIndex(rootIndex),
      _existingNodePtrMap({
        {startingFrom, {_rootIndex}}
      })
  {}

  template<typename Graph>
  void discover_vertex(const AtomIndexType& u, const Graph&) {
    assert(_depthMap.count(u));
  }

  /*!
   * This is called after all edges pertaining to this vertex are called, and
   * always before descending to a new search depth. This method places all
   * backwards edges that were remembered in non-tree-edge calls on all
   * available nodes of the same parent.
   */
  template<typename Graph>
  void finish_vertex(const AtomIndexType&, const Graph&) {
    auto& tree = _treeBaseRef._tree;

    TemplateMagic::inplaceRemoveIf(
      _edgesToAddBeforeDescent,
      [&](const auto& edgePair) {
        bool canBeAdded = (
          _existingNodePtrMap.count(edgePair.first) > 0
          && _existingNodePtrMap.count(edgePair.second) > 0
        );

        // For every node with the target index
        if(canBeAdded) {
          for(const auto& parentNodeIndex : _existingNodePtrMap.at(edgePair.second)) {
            /* if the depth of the target node in the generated tree is
             * less or equal to the depth of the source vertex in the graph,
             * then add this new node as child
             */
            if(
              _treeBaseRef._depthOfNode(parentNodeIndex) <= _depthMap[edgePair.first]
              && !tree[parentNodeIndex].isDuplicate
            ) {
              // Create a new child node and add it to the parent
              auto newNodeIndex = boost::add_vertex(tree);
              tree[newNodeIndex].molIndex = edgePair.first;

              boost::add_edge(parentNodeIndex, newNodeIndex, tree);

              tree[newNodeIndex].isDuplicate = _treeBaseRef._molIndexExistsInBranch(
                edgePair.first, 
                parentNodeIndex
              );

              _treeBaseRef._addBondOrderDuplicates(parentNodeIndex, newNodeIndex);

              // Add it to existingNodePtrMap
              if(_existingNodePtrMap.count(edgePair.first) == 0) {
                _existingNodePtrMap[edgePair.first] = {newNodeIndex};
              } else {
                _existingNodePtrMap[edgePair.first].push_back(newNodeIndex);
              }
            }
          }
        }

        return canBeAdded;
      }
    );
  }

  /*!
   * We're not paricularly interested in tree edges aside from ensuring that
   * our depth map stays correct.  Why do we do (basically) nothing on tree
   * edges? Because we want to add edges ONLY after we have completely
   * discovered the sphere that contains the source edges here
   */
  template<typename Graph>
  void tree_edge(const EdgeIndexType& e, const Graph& g) {
    /* 
     */
    // Assign the depth map
    _depthMap[boost::target(e, g)] = _depthMap[boost::source(e, g)] + 1;
  }

  /*! 
   * This is called on all cycle closing edges and backwards-pointing edges,
   * i.e. to the previous search depth. 
   * 
   * We act on non-tree edges because precisely then are the
   * backward-pointing edges shown, up into the previous shell of the BFS
   * iteration, which is now fully expanded and stored in the depth map.
   */
  template<typename Graph>
  void non_tree_edge(const EdgeIndexType& e, const Graph& g) {
    auto& tree = _treeBaseRef._tree;

    auto source = boost::source(e, g), // where we are now
         target = boost::target(e, g); // maybe points to lower depth vertex

    /* Add the source as child of target, provided target is in the tree and
     * is at a lower depth than the source 
     */
    if( 
      /* Is this truly a backward-pointing edge, i.e. it points from the
       * current exploration shell of BFS further inwards?
       */
      _existingNodePtrMap.count(target)
      && _depthMap[target] < _depthMap[source]
    ) {
      // For every node with the target index
      for(const auto& parentNodeIndex : _existingNodePtrMap.at(target)) {
        if(!tree[parentNodeIndex].isDuplicate) {
          auto newNodeIndex = boost::add_vertex(tree);
          tree[newNodeIndex].molIndex = source;
          boost::add_edge(parentNodeIndex, newNodeIndex, tree);

          tree[newNodeIndex].isDuplicate = _treeBaseRef._molIndexExistsInBranch(
            source, 
            parentNodeIndex
          );

          _treeBaseRef._addBondOrderDuplicates(parentNodeIndex, newNodeIndex);

          // Add it to existingNodePtrMap
          if(_existingNodePtrMap.count(source) == 0) {
            _existingNodePtrMap[source] = {newNodeIndex};
          } else {
            _existingNodePtrMap[source].push_back(newNodeIndex);
          }
        }
      }
    } else if(_depthMap[target] == _depthMap[source]) {
      /* If the edge is on the same shell, then it must cycle-closing, so we
       * have to wait until all vertices in this depth are discovered and only
       * prior to descent add them.
       */
      _edgesToAddBeforeDescent.emplace_back(source, target);
    }
  }
};

RankingTree::RankingTree(
  const Molecule& molecule,
  const AtomIndexType& atomToRank
) : _moleculeRef(molecule) {
  /* Part 1: Make a complete acyclic graph */

  auto rootIndex = boost::add_vertex(_tree);
  _tree[rootIndex].molIndex = atomToRank;
  _tree[rootIndex].isDuplicate = false;

  using ColorMapBase = std::map<AtomIndexType, boost::default_color_type>;

  ColorMapBase colorMap;
  boost::associative_property_map<ColorMapBase> propColorMap(colorMap);
  boost::queue<GraphType::vertex_descriptor> Q;

  RankingTreeGenerator visitor(
    *this,
    atomToRank,
    rootIndex
  );

  try {
    boost::breadth_first_visit(
      // The graph to operate on
      molecule.getGraph(),
      // The vertex to start with
      atomToRank,
      // A queue object to store vertex_descriptors
      Q,
      // The visitor to use
      visitor,
      // A map to store color (state)
      propColorMap
    );
  } catch(RankingTreeGenerator::EarlyExit& e) {}

  // Now that the graph is acyclic, finish the tree by DFS at every leaf
  TreeGraphType::vertex_iterator iter, end;
  std::tie(iter, end) = boost::vertices(_tree);

  // Skip the root vertex, that one does not need DFS-ing
  ++iter;

  while(iter != end) {
    if(!_tree[*iter].isDuplicate) {
      /* We can insert vertices and edges while iterating through, since
       * these do not invalidate vertex iterators
       * (boost.org/doc/libs/1_65_1/libs/graph/doc/adjacency_list.html)
       */
      _DFSExpand(*iter, _molIndicesInBranch(*iter));
    }

    ++iter;
  }
}

class RankingTree::GraphvizWriter {
private:
  // Closures
  const RankingTree& _baseRef;
  const std::string _title;
  const std::set<TreeVertexIndex> _squareVertices;
  const std::set<TreeVertexIndex> _colorVertices;
  const std::set<TreeEdgeIndex> _colorEdges;

public:
  GraphvizWriter(
    const RankingTree& baseTree,
    const std::string& title = "",
    const std::set<TreeVertexIndex>& squareVertices = {},
    const std::set<TreeVertexIndex>& colorVertices = {},
    const std::set<TreeEdgeIndex>& colorEdges = {}
  ) : _baseRef(baseTree),
      _title(title),
      _squareVertices(squareVertices),
      _colorVertices(colorVertices),
      _colorEdges(colorEdges)
  {}

  void operator() (std::ostream& os) const {
    os << R"(  graph [fontname = "Arial"];)" << "\n"
      << R"(  node [fontname = "Arial", shape = "circle", style = "filled"];)" << "\n"
      << R"(  edge [fontname = "Arial"];)" << "\n"
      << R"(  labelloc="t"; label=")" << _title << R"(")" << ";\n";
  }

  void operator() (std::ostream& os, const TreeVertexIndex& vertexIndex) const {
    auto symbolString = Delib::ElementInfo::symbol(
      _baseRef._moleculeRef.getElementType(
        _baseRef._tree[vertexIndex].molIndex
      )
    );

    bool hasStereocenter = _baseRef._tree[vertexIndex].stereocenterOption.operator bool();

    os << "[" 
      << R"(label=")" << vertexIndex << "-" << symbolString 
      << _baseRef._tree[vertexIndex].molIndex << R"(")";

    // Node background coloring
    if(_colorVertices.count(vertexIndex) > 0) {
      os << R"(, fillcolor="tomato")";
    } else if(hasStereocenter) {
      os << R"(, fillcolor="steelblue")";
    } else if(MolGraphWriter::elementBGColorMap.count(symbolString) != 0u) {
      os << R"(, fillcolor=")" 
        << MolGraphWriter::elementBGColorMap.at(symbolString) << R"(")";
    } 

    // Font coloring
    if(_colorVertices.count(vertexIndex) > 0) {
      os << R"(, fontcolor="white")";
    } else if(MolGraphWriter::elementTextColorMap.count(symbolString) != 0u) {
      os << R"(, fontcolor=")" 
        << MolGraphWriter::elementTextColorMap.at(symbolString) << R"(")";
    } else if(hasStereocenter) {
      os << R"(, fontcolor="white")";
    }

    // Shape alterations
    if(_squareVertices.count(vertexIndex) > 0) {
      os << R"(, shape="square")";
    } else if(_baseRef._tree[vertexIndex].isDuplicate) {
      // Duplicate atoms get double-circle shape
      os << R"(, shape="doublecircle")";
    } else if(hasStereocenter) {
      os << R"(, shape="diamond")";
    }

    // Tooltip
    if(hasStereocenter) {
      os << R"(, tooltip=")" 
        << _baseRef._tree[vertexIndex].stereocenterOption.value().info()
        << R"(")";
    }

    // Make hydrogens smaller
    if(symbolString == "H") {
      os << ", fontsize=10, width=.6, fixedsize=true";
    }

    os << "]";
  }

  void operator() (std::ostream& os, const TreeEdgeIndex& edgeIndex) const {
    auto hasStereocenter = _baseRef._tree[edgeIndex].stereocenterOption.operator bool();

    os << "[";

    // Coloring
    if(_colorEdges.count(edgeIndex) > 0) {
      os << R"(color="tomato")";
    } else if(hasStereocenter) {
      os << R"(color="steelblue")";
    }
    
    // Edge width
    if(_colorEdges.count(edgeIndex) > 0 || hasStereocenter) {
      os << R"(, penwidth="2")";
    }

    // Tooltip
    if(hasStereocenter) {
      os << R"(, tooltip=")"
        << _baseRef._tree[edgeIndex].stereocenterOption.value().info()
        << R"(")";
    }

    os << "]";
  }
};

/* Sequence rule classes */
struct RankingTree::SequenceRuleOneVertexComparator {
private:
  const RankingTree& _base;

public:
  SequenceRuleOneVertexComparator(const RankingTree& base) : _base(base) {}

  /* NOTE: Since the multisets this is used in need a DESC ordering,
   * The arguments names have been swapped, effectively inverting the sorting
   */
  bool operator () (const TreeVertexIndex& b, const TreeVertexIndex& a) const {
    /* Cases:
     * - Neither is duplicate -> Compare Zs
     * - One is duplicate -> Non-duplicate is bigger
     * - Both are duplicate -> Compare distance to root
     */
    const bool& aDuplicate = _base._tree[a].isDuplicate;
    const bool& bDuplicate = _base._tree[b].isDuplicate;

    if(!aDuplicate && !bDuplicate) {
      // Casting elementType to unsigned basically gives Z
      return static_cast<unsigned>(
        _base._moleculeRef.getElementType(_base._tree[a].molIndex)
      ) < static_cast<unsigned>(
        _base._moleculeRef.getElementType(_base._tree[b].molIndex)
      );
    } 
    
    if(aDuplicate && bDuplicate) {
      // That vertex is smaller which is further from root (invert comparison)
      return _base._duplicateDepth(a) > _base._duplicateDepth(b);
    }

    if(aDuplicate && !bDuplicate) {
      return true;
    }

    // The case below is not needed, it falls together with equality
    /*if(!aDuplicate && bDuplicate) {
      return false;
    }*/

    return false;
  }
};

struct RankingTree::SequenceRuleThreeEdgeComparator {
private:
  const RankingTree& _base;

public:
  SequenceRuleThreeEdgeComparator(const RankingTree& base) : _base(base) {}

  bool operator () (const TreeEdgeIndex& a, const TreeEdgeIndex& b) const {
    /* Cases:
     * - Neither is instantiated -> neither can be smaller
     * - One is instantiated -> that one is bigger
     * - Both are instantiated -> Compare sequentially:
     *   - Number of possible assignments (non-stereogenic is smaller)
     *   - If assigned or not (unassigned is smaller)
     *   - Assignment value (Z > E)
     *
     * NOTE: Since the multisets this is used in need a DESC ordering,
     * all comparisons below are reversed.
     */

    auto EZOptionalA = _base._tree[a].stereocenterOption;
    auto EZOptionalB = _base._tree[b].stereocenterOption;

    if(!EZOptionalA && !EZOptionalB) {
      /* This does not invalidate the program, it just means that all
       * uninstantiated EZStereocenters are equal, and no relative ordering
       * is provided.
       */
      return false;
    }

    if(EZOptionalA && !EZOptionalB) {
      return true;
    }

    if(!EZOptionalA && EZOptionalB) {
      return false;
    }

    // Now we know both actually have a stereocenterPtr
    const auto& EZStereocenterA = EZOptionalA.value();
    const auto& EZStereocenterB = EZOptionalB.value();

    // Reverse everything below for descending sorting
    return TemplateMagic::componentSmaller(
      EZStereocenterB.numAssignments(),
      EZStereocenterA.numAssignments()
    ).value_or(
      // Mixed optional comparison (includes comparison of assignment value if assigned)
      EZStereocenterB.assigned() < EZStereocenterA.assigned()
    );
  }
};

struct RankingTree::SequenceRuleFourVariantComparator {
public:
  using CompareType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

  class VariantComparisonVisitor : boost::static_visitor<bool> {
  private:
    const TreeGraphType& _treeRef;

  public:
    explicit VariantComparisonVisitor(
      const SequenceRuleFourVariantComparator& comparatorBase
    ) : _treeRef(comparatorBase._base._tree) {}

    /* Surprisingly, the code for edge and vertex indices is completely
     * identical, so we can abstract over the types.
     *
     * NOTE: Since we desire the inverted ordering sequence, this comparison
     * is implemented as normal, but swapping the parameters a and b.
     */
    template<typename T>
    bool operator() (const T& b, const T& a) const {
      const auto& aOption = _treeRef[a].stereocenterOption;
      const auto& bOption = _treeRef[b].stereocenterOption;

      // Uninstantiated stereocenters always compare false
      if(!aOption && !bOption) {
        return false;
      }

      // Instantiated is less than uninstantiated
      if(aOption && !bOption) {
        return false;
      }

      if(!aOption && bOption) {
        return true;
      }

      // Now we know both actually have an instantiated stereocenter
      const auto& StereocenterA = aOption.value();
      const auto& StereocenterB = bOption.value();

      /* We only want to know whether they differ in stereogenicity, i.e.
       * whether one is stereogenic while the other isn't. In effect, this
       * means comparing whether they have more than one assignment:
       *
       * Cases:
       * - Neither is stereogenic -> false
       * - Both are stereogenic -> false (neither bigger than other)
       * - A isn't stereogenic, B is
       *   false < true == 0 < 1 == true
       *   (So A < B is true, meaning stereogenic < non-stereogenic, leading
       *   to the desired DESC ordering)
       *
       * This is valid for both CN and EZ types of stereocenters
       */
      return (
        (StereocenterA.numAssignments() > 1) 
        < (StereocenterB.numAssignments() > 1)
      );
    }

    // For different types
    template<typename T, typename U>
    bool operator() (const T&, const U&) const {
      return false;
    }
  };

private:
  const RankingTree& _base;
  VariantComparisonVisitor _variantComparator;

public:
  SequenceRuleFourVariantComparator(const RankingTree& base) 
    : _base(base),
      _variantComparator {*this}
  {}

  bool operator () (const CompareType& a, const CompareType& b) const {
    return boost::apply_visitor(_variantComparator, a, b);
  }
};

struct RankingTree::SequenceRuleFiveVariantComparator {
public:
  using CompareType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

  class VariantComparisonVisitor : boost::static_visitor<bool> {
  private:
    const TreeGraphType& _treeRef;

  public:
    explicit VariantComparisonVisitor(
      const SequenceRuleFiveVariantComparator& base
    ) : _treeRef(base._base._tree) {}

    /* The code for edge and vertex indices is completely identical, so we
     * can abstract over the types.
     *
     * NOTE: Since we desire the inverted ordering sequence, this comparison
     * is implemented as normal, but swapping the parameters a and b.
     */
    template<typename T>
    bool operator() (const T& b, const T& a) const {
      const auto& aOption = _treeRef[a].stereocenterOption;
      const auto& bOption = _treeRef[b].stereocenterOption;

      // Uninstantiated stereocenters always compare false
      if(!aOption && !bOption) {
        return false;
      }

      // Instantiated is less than uninstantiated
      if(aOption && !bOption) {
        return false;
      }

      if(!aOption && bOption) {
        return true;
      }

      // Now we know both actually have an instantiated stereocenter
      const auto& StereocenterA = aOption.value();
      const auto& StereocenterB = bOption.value();

      // Need to compare using the number of assignments first
      if(StereocenterA.numAssignments() < StereocenterB.numAssignments()) {
        return true;
      }

      if(StereocenterB.numAssignments() < StereocenterA.numAssignments()) {
        return false;
      }

      /* In our context, sequence rule 5 directly compares the assignment
       * of the assigned stereocenters.
       */
      return StereocenterA.assigned() < StereocenterB.assigned();
    }

    // For different types
    template<typename T, typename U>
    bool operator() (const T&, const U&) const {
      return false;
    }
  };

private:
  const RankingTree& _base;
  VariantComparisonVisitor _variantComparator;

public:
  SequenceRuleFiveVariantComparator(const RankingTree& base) 
    : _base(base),
      _variantComparator {*this}
  {}

  bool operator () (const CompareType& a, const CompareType& b) const {
    return boost::apply_visitor(_variantComparator, a, b);
  }
};

/* Variant visitor helper classes */
class RankingTree::VariantHasInstantiatedStereocenter : boost::static_visitor<bool> {
private:
  const RankingTree& _baseRef;

public:
  explicit VariantHasInstantiatedStereocenter(
    const RankingTree& base
  ) : _baseRef(base) {}

  //! Check if the variant is an instantiated stereocenter
  template<typename T>
  bool operator() (const T& a) const {
    const auto& stereocenterOption = _baseRef._tree[a].stereocenterOption;
    return stereocenterOption.operator bool();
  }
};

class RankingTree::VariantDepth : boost::static_visitor<unsigned> {
private:
  const RankingTree& _baseRef;

public:
  explicit VariantDepth(
    const RankingTree& base
  ) : _baseRef(base) {}

  //! Returns the depth of the vertex or edge
  template<typename T>
  unsigned operator() (const T& a) const {
    return _baseRef._mixedDepth(a);
  }
};

class RankingTree::VariantSourceNode : boost::static_visitor<TreeVertexIndex> {
private:
  const RankingTree& _baseRef;

public:
  explicit VariantSourceNode(
    const RankingTree& base
  ) : _baseRef(base) {}

  //! Returns the depth of the vertex or edge
  TreeVertexIndex operator() (const TreeVertexIndex& a) const {
    return a;
  }

  TreeVertexIndex operator() (const TreeEdgeIndex& a) const {
    return boost::source(a, _baseRef._tree);
  }
};

class RankingTree::VariantStereocenterStringRepresentation : boost::static_visitor<std::string> {
private:
  const RankingTree& _baseRef;

public:
  explicit VariantStereocenterStringRepresentation(
    const RankingTree& base
  ) : _baseRef(base) {}

  //! Returns a string representation of the *type* of the stereocenter
  template<typename T>
  std::string operator() (const T& a) const {
    const auto& aOption = _baseRef._tree[a].stereocenterOption;
    if(aOption) {
      return aOption.value().rankInfo();
    }

    return "";
  }
};

class RankingTree::VariantLikePair : boost::static_visitor<bool> {
private:
  const RankingTree& _baseRef;

public:
  explicit VariantLikePair(const RankingTree& base) : _baseRef(base) {}

  //! Returns a string representation of the *type* of the stereocenter
  template<typename T, typename U>
  bool operator() (const T& a, const U& b) const {
    const auto& aOption = _baseRef._tree[a].stereocenterOption;
    const auto& bOption = _baseRef._tree[b].stereocenterOption;

    return (
      aOption 
      && bOption
      && aOption.value().numAssignments() == bOption.value().numAssignments()
      && aOption.value().assigned() == bOption.value().assigned()
    );
  }
};


/*!
 * This function ranks the direct substituents of the atom it is instantiated
 * upon by the sequential application of the 2013 IUPAC Blue Book sequence
 * rules. They are somewhat adapted since the priority of asymmetric centers
 * of higher (and also lower) symmetry must also be considered because
 * transition metal chemistry is also included in this library.
 *
 * It returns a sorted vector of vectors, in which every sub-vector
 * represents a set of equal-priority substituents. The sorting is ascending,
 * meaning from lowest priority to highest priority.
 */
std::vector<
  std::vector<AtomIndexType>
> RankingTree::rank(const std::set<AtomIndexType>& excludeIndices) {

  /* Find out how direct children of the root node are ordered */

  /* Represent overall relationship of root's direct children vertex indices
   * to be ranked against one another
   */
  // Get all direct children of the root node, which are specifically non-duplicate
  std::set<TreeVertexIndex> rootChildren;
  {
    TreeGraphType::out_edge_iterator iter, end;
    std::tie(iter, end) = boost::out_edges(0, _tree);

    while(iter != end) {
      auto childIndex = boost::target(*iter, _tree);

      /* Limit the ranking of root's children to non-duplicate and
       * non-excluded nodes
       */
      if(
        !_tree[childIndex].isDuplicate 
        && excludeIndices.count(_tree[childIndex].molIndex) == 0
      ) {
        rootChildren.insert(childIndex);
      }

      ++iter;
    }
  }

  /* Create a helper object with which we can gradually discover the ordering
   * relationships
   */
  auto orderingHelper = OrderDiscoveryHelper<TreeVertexIndex>(rootChildren);


#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "Ranking substituents of atom index " 
    << _tree[0].molIndex
    << ": {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

  /* Prior to returning, we must transform the ranked tree vertex indices
   * back to molecule indices. This pattern is needed after the full
   * application of every sequence rule, so it is declared here as a lambda.
   */
  auto fetchResult = [&]() -> std::vector<
    std::vector<AtomIndexType>
  > {
    return TemplateMagic::map(
      orderingHelper.getSets(),
      [&](const auto& set) -> std::vector<AtomIndexType> {
        return TemplateMagic::map(
          set,
          [&](const auto& treeVertex) -> AtomIndexType {
            return _tree[treeVertex].molIndex;
          }
        );
      }
    );
  };
#ifndef NDEBUG 
  auto visitedVertices = [](
    const std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    >& seeds,
    const std::vector<
      std::vector<TreeVertexIndex>
    >& undecidedSets
  ) -> std::set<TreeVertexIndex> {
    std::set<TreeVertexIndex> visited;

    for(const auto& undecidedSet : undecidedSets) {
      for(const auto& undecidedBranch : undecidedSet) {
        for(const auto& seed : seeds.at(undecidedBranch)) {
          visited.insert(seed);
        }
      }
    }

    return visited;
  };
#endif

  /* A much-needed variable which was often re-declared within the local
   * scopes of every sequence rule. This allows reuse.
   */
  auto undecidedSets = orderingHelper.getUndecidedSets();

  /* Evaluate by sequence rule 1
   * - Higher atomic number precedes lower atomic number
   * - A duplicate atom node whose corresponding non-duplicate atom node is
   *   root or is closer to the root preceds a duplicate node whose
   *   corresponding atom node is further from the root
   */
  { // New scope to avoid namespace pollution
    /* For each branch to compare, keep a multiset of all the encountered
     * indices in the branch. The multiset keeps a descending order of the
     * indices according to sequence rule 1 and can be lexicographically
     * compared against one another using the SequenceRuleOneVertexComparator
     *
     * NOTE: Do not use operator [] to access members in comparisonSets. This
     * old form of access and assignable proxy object forces a
     * default-instantiation of the Comparator, which
     * SequenceRuleOneVertexComparator cannot do. Always use .at().
     */
    std::map<
      TreeVertexIndex,
      std::multiset<TreeVertexIndex, SequenceRuleOneVertexComparator>
    > comparisonSets;

    for(const auto& childIndex : rootChildren) {
      comparisonSets.emplace(
        childIndex, // KeyType=TreeVertexIndex ctor
        *this // MappedType=multiset ctor requires a reference to *this
      );

      comparisonSets.at(childIndex).insert(childIndex);
    }

    /* For each branch to compare, keep a set of seed indices to follow in an
     * iteration. These are expanded in each iteration as long as the branch
     * they contain remains relevant (i.e. it's relation to the other branches
     * is still unclear): They themselves are placed in the comparisonSet, 
     * while their respective children are the new seeds for the next
     * iteration.
     */
    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    for(const auto& childIndex : rootChildren) {
      seeds.emplace(
        childIndex,
        _children(childIndex)
      );
    }

    // Perform the first comparison, across ALL root children
    TemplateMagic::forAllPairs(
      rootChildren,
      [&](const TreeVertexIndex& a, const TreeVertexIndex& b) {
        if(_multisetCompare(comparisonSets.at(a), comparisonSets.at(b))) {
          orderingHelper.addLessThanRelationship(a, b);
        } else if(_multisetCompare(comparisonSets.at(b), comparisonSets.at(a))) {
          orderingHelper.addLessThanRelationship(b, a);
        }
      }
    );

    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 1", {0}, visitedVertices(seeds, undecidedSets)),
        _makeGraph("Sequence rule 1 multisets", 0, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    /* Checks if there are seeds whose expansion could be relevant for
     * undecided indices
     */
    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a full BFS Step on all undecided set seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {
            // Add the seed to the indices used for comparison
            comparisonSets.at(undecidedTreeIndex).insert(seed);

            // Add its children to the new seeds
            for(const auto& childIndex : _children(seed)) {
              newSeeds.insert(childIndex);
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 1", {0}, visitedVertices(seeds, undecidedSets)),
          _makeGraph("Sequence rule 1 multisets", 0, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }
  }

#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "Sets post sequence rule 1: {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

  // Was sequence rule 1 enough?
  if(orderingHelper.isTotallyOrdered()) {
    return fetchResult();
  }

  /* Sequence rule 2
   * - A node with higher atomic mass precedes ones with lower atomic mass
   *
   * NOTE: We skip this sequence rule for two reasons
   * - There is currently no way to represent isotopes in the library
   * - This rule may come to bear when disconnecting heteroaromatic cycles,
   *   where all aromatic atoms are given a duplicate atom with a mass equal
   *   to the average of what it would have if the double bonds were located
   *   at each of the possible positions. Although differentiating these
   *   cases would be good, there is currently no code to detect aromaticity
   *   nor permute all possible arrangements of double bonds in their
   *   subgraphs.
   */


  /* Sequence rule 3: double bond configurations
   * - Z > E > unassigned > non-stereogenic 
   *   (seqCis > seqTrans > non-stereogenic)
   */
  /* Before we can compare using sequence rule 3, we have to instantiate all
   * auxiliary descriptors, meaning we have to instantiate EZStereocenters in
   * all double bond positions and CNStereocenters in all candidate
   * positions. 
   *
   * In both cases, we have to rank non-duplicate adjacent vertices according
   * to the all of the same sequence rules, with the difference that graph
   * traversal isn't straightforward top-down, but proceeds in all unexplored
   * directions. 
   *
   * We want to proceed from the outer parts of the tree inwards, so that
   * branches that are forced to consider later sequence rules have their
   * required stereocenters instantiated already.
   */
  // Cleanly divide the tree into depths by its edges
  std::vector<
    std::set<TreeEdgeIndex>
  > byDepth;

  { // Populate byDepth from active branches only
    std::set<TreeEdgeIndex> rootEdges;

    for(const auto& undecidedSet : undecidedSets) {
      for(const auto& undecidedBranch : undecidedSet) {
        rootEdges.insert(
          boost::edge(0, undecidedBranch, _tree).first
        );
      }
    }

    byDepth.emplace_back(
      std::move(rootEdges)
    );
  }

  std::set<TreeEdgeIndex> nextChildren;

  for(const auto& edge : byDepth.back()) {
    auto outIterPair = boost::out_edges(
      boost::target(edge, _tree), _tree
    );

    while(outIterPair.first != outIterPair.second) {
      nextChildren.insert(*outIterPair.first);

      ++outIterPair.first;
    }
  }

  while(nextChildren.size() != 0) {
    byDepth.emplace_back(
      std::move(nextChildren)
    );

    nextChildren.clear();

    for(const auto& edge : byDepth.back()) {
      auto outIterPair = boost::out_edges(
        boost::target(edge, _tree), _tree
      );

      while(outIterPair.first != outIterPair.second) {
        nextChildren.insert(*outIterPair.first);

        ++outIterPair.first;
      }
    }
  }

#ifndef NDEBUG
  auto collectHandledEdges = [&](
    const typename decltype(byDepth)::reverse_iterator& rIter
  ) -> std::set<TreeEdgeIndex> {
    std::set<TreeEdgeIndex> edgeIndices;

    for(auto it = byDepth.rbegin(); it != rIter; ++it) {
      const auto& currentEdges = *it;

      for(const auto& edge : currentEdges) {
        edgeIndices.insert(edge);
      }
    }

    return edgeIndices;
  };
#endif

  // Remember if you found something to instantiate or not
  bool foundStereocenters = false;

  // Process the tree, from the bottom up
  for(auto it = byDepth.rbegin(); it != byDepth.rend(); ++it) {
    const auto& currentEdges = *it;

    for(const auto& edge: currentEdges) {
      auto sourceIndex = boost::source(edge, _tree);

      if(_doubleBondEdges.count(edge)) {
        // Instantiate an EZStereocenter here!
        auto targetIndex = boost::target(edge, _tree);

        RankingInformation sourceRanking, targetRanking;

        // Get ranking for the edge substituents
        sourceRanking.sortedSubstituents = _omniDirectionalRank(
          sourceIndex,
          {targetIndex}
        );
        targetRanking.sortedSubstituents = _omniDirectionalRank(
          targetIndex,
          {sourceIndex}
        );

        /* NOTE: There is no need to collect linking information since we
         * are in an acyclic digraph, and any cycles in the original molecule
         * are no longer present.
         */

        _tree[edge].stereocenterOption = Stereocenters::EZStereocenter {
          sourceIndex,
          sourceRanking,
          targetIndex,
          targetRanking
        };

        // Mark that we instantiated something
        foundStereocenters = true;

#ifndef NDEBUG
        if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
          _writeGraphvizFiles({
            _adaptMolGraph(_moleculeRef.dumpGraphviz()),
            dumpGraphviz("Sequence rule 3 preparation", {0}, {}, collectHandledEdges(it))
          });
        }
#endif
      }

      if(
        !_tree[edge].stereocenterOption // No EZStereocenter on this edge
        && !_tree[sourceIndex].stereocenterOption // No CNStereocenter
        && _nonDuplicateDegree(sourceIndex) >= 3 // Min. degree for chirality
        && sourceIndex != 0 // not root, no CNStereocenters needed there!
      ) {
        /* TODO reduce the effort here by counting the number of adjacent
         * (terminal!) hydrogens and comparing with a table of symmetries vs
         * #hydrogens 
         * -> max # assignments
         *
         * In case only one assignment is possible, there is no reason to rank
         * the substituents
         */
        // Instantiate a CNStereocenter here!
        RankingInformation centerRanking;
        centerRanking.sortedSubstituents = _omniDirectionalRank(
          sourceIndex,
          {}
        );

        _tree[sourceIndex].stereocenterOption = Stereocenters::CNStereocenter {
          _moleculeRef.determineLocalGeometry(
            _tree[sourceIndex].molIndex
          ),
          sourceIndex,
          centerRanking
        };

        // Mark that we instantiated something
        foundStereocenters = true;
        
#ifndef NDEBUG
        if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
          _writeGraphvizFiles({
            _adaptMolGraph(_moleculeRef.dumpGraphviz()),
            dumpGraphviz("Sequence rule 3 preparation", {0}, {}, collectHandledEdges(it))
          });
        }
#endif
      }

    }
  }

  // Was anything instantiated? If not, we can skip rules 3 through 5.
  if(!foundStereocenters) {
    return fetchResult();
  }

  // Now, we can apply sequence rule 3
  { // Sequence rule 3 local scope
    // Edge BFS outwards from each
    std::map<
      TreeVertexIndex,
      std::multiset<TreeEdgeIndex, SequenceRuleThreeEdgeComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    undecidedSets = orderingHelper.getUndecidedSets();

    for(const auto& undecidedSet : undecidedSets) {
      for(const auto& undecidedBranch : undecidedSet) {
        comparisonSets.emplace(
          undecidedBranch,
          *this
        );

        auto startingEdge = boost::edge(0, undecidedBranch, _tree);
        assert(startingEdge.second);

        comparisonSets.at(undecidedBranch).insert(
          startingEdge.first
        );

        seeds[undecidedBranch] = {undecidedBranch};
      }
    }

    /* Perform first comparison (undecided sets up-to-date from previous
     * sequence rule)
     */
    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    // Re-calculate the undecided sets
    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 3", {0}, visitedVertices(seeds, undecidedSets)),
        _makeGraph("Sequence rule 3 multisets", 0, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a BFS step on all seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {
            // Add all adjacent out-edges on seed to the comparison multiset
            for(
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdgeIndex = *outIterPair.first;

              comparisonSets.at(undecidedTreeIndex).insert(outEdgeIndex);
              
              // Add the target node as a seed only if it's non-terminal
              auto targetIndex = boost::target(outEdgeIndex, _tree);
              if(boost::out_degree(targetIndex, _tree) > 0) {
                newSeeds.insert(targetIndex);
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 3", {0}, visitedVertices(seeds, undecidedSets)),
          _makeGraph("Sequence rule 3 multisets", 0, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }
  } // End sequence rule 3 local scope

#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "Sets post sequence rule 3: {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

  // Was sequence rule 3 enough?
  if(orderingHelper.isTotallyOrdered()) {
    return fetchResult();
  }

  /* Sequence rule 4: Other configurations (RS, MP, EZ)
   * - stereogenic > pseudostereogenic > non-stereogenic
   * - like pairs of descriptors over unlike pairs of descriptors
   *   set 1: {R, M, Z}, set 2: {S, P, E}, like within a set, unlike crossover
   *
   *   Procedure:
   *
   *   - Create a set of all stereogenic groups within the branches to rank
   *   - Establish relative rank by application of sequence rules. 
   *
   *     Can probably use _omniDirectionalRank recursively at the first
   *     junction between competing stereocenters on the respective branch
   *     indices to find out which precedes which, provided both are at the
   *     same depth from root.
   *
   *     (WARNING: This may be tricky!)
   *
   *   - Choose a representative stereodescriptor for each branch according to
   *     the following rules:
   *
   *     - The solitarily highest ranked stereogenic group
   *     - The stereodescriptor occurring more often than all others in 
   *       the set of stereodescriptors
   *     - All most-often-occurring stereodescriptors
   *
   *   - Precedence:
   *     - The branch with fewer representative stereodescriptors has priority
   *     - Go through both lists of ranked stereocenters simultaneously,
   *       proceeding from highest rank to lowest, pairing each occurring
   *       descriptor with the reference descriptor.
   *
   *       like pairs precede unlike pairs of stereodescriptors
   *
   *
   * - (Pseudoasymmetries) r precedes s and m precedes p
   */

  { // Sequence rule 4 local scope
    /* First part of rule: A regular omni-directional BFS, inserting both
     * edges and node indices into the multiset! Perhaps avoiding the
     * insertion of uninstantiated edges and vertices entirely?
     */

    using VariantType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

    std::map<
      TreeVertexIndex,
      std::multiset<VariantType, SequenceRuleFourVariantComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    // Initialize the BFS state
    for(const auto& undecidedSet : orderingHelper.getUndecidedSets()) {
      for(const auto& undecidedBranch : undecidedSet) {
        comparisonSets.emplace(
          undecidedBranch,
          *this
        );

        comparisonSets.at(undecidedBranch).emplace(undecidedBranch);
        comparisonSets.at(undecidedBranch).emplace(
          boost::edge(0, undecidedBranch, _tree).first
        );

        seeds[undecidedBranch] = {undecidedBranch};
      }
    }

    /* First comparison (undecided sets are up-to-date from previous sequence
     * rule)
     */
    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    // Recalculate undecided sets
    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 4A", {0}, visitedVertices(seeds, undecidedSets)),
        _makeGraph("Sequence rule 4A multisets", 0, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a full BFS Step on all undecided set seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {
            for( // Out-edges
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdge = *outIterPair.first;

              auto edgeTarget = boost::target(outEdge, _tree);

              if(
                comparisonSets.at(undecidedTreeIndex).count(edgeTarget) == 0
                && !_tree[edgeTarget].isDuplicate
              ) {
                comparisonSets.at(undecidedTreeIndex).emplace(edgeTarget);
                comparisonSets.at(undecidedTreeIndex).emplace(outEdge);

                newSeeds.insert(edgeTarget);
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 4A", {0}, visitedVertices(seeds, undecidedSets)),
          _makeGraph("Sequence rule 4A multisets", 0, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }

#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "Sets post sequence rule 4A: {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

    // Is Sequence Rule 4, part A enough?
    if(orderingHelper.isTotallyOrdered()) {
      return fetchResult();
    }

    /* Second part: Choosing representational stereodescriptors for all
     * branches and pairing in sequence of priority for branches.
     */
    std::map<
      TreeVertexIndex,
      std::set<VariantType>
    > stereocenterMap;

    VariantHasInstantiatedStereocenter isInstantiatedChecker {*this};

    // Copy all those variants from part A that are actually instantiated
    for(const auto& undecidedSet : orderingHelper.getUndecidedSets()) {
      for(const auto& undecidedBranch : undecidedSet) {
        stereocenterMap[undecidedBranch] = {};

        for(const auto& variantValue : comparisonSets.at(undecidedBranch)) {
          if(boost::apply_visitor(isInstantiatedChecker, variantValue)) {
            stereocenterMap.at(undecidedBranch).insert(variantValue);
          }
        }
      }
    }

    std::map<
      TreeVertexIndex,
      OrderDiscoveryHelper<VariantType>
    > relativeOrders;

    VariantDepth depthFetcher {*this};
    VariantSourceNode sourceNodeFetcher {*this};
    VariantStereocenterStringRepresentation stringRepFetcher {*this};

    std::map<
      TreeVertexIndex,
      std::set<VariantType>
    > representativeStereodescriptors;

    for(const auto& mapIterPair : stereocenterMap) {
      const auto& branchIndex = mapIterPair.first;
      const auto& variantSet = mapIterPair.second;

      // Instantiate the order discovery helpers
      relativeOrders.emplace(
        branchIndex,
        variantSet
      );

      // Compare variants in a branch based on mixed depth
      TemplateMagic::forAllPairs(
        variantSet,
        [&](const auto& a, const auto& b) {
          auto aDepth = boost::apply_visitor(depthFetcher, a);
          auto bDepth = boost::apply_visitor(depthFetcher, b);

          if(aDepth < bDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
          } else if(bDepth < aDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
          }
        }
      );

      /* Try to resolve undecided sets via ranking downwards at junction of
       * pairs
       */
      for(
        const auto& undecidedSet : 
        relativeOrders.at(branchIndex).getUndecidedSets()
      ) {
        TemplateMagic::forAllPairs(
          undecidedSet,
          [&](const auto& a, const auto& b) {
            auto junctionInfo = _junction(
              boost::apply_visitor(sourceNodeFetcher, a),
              boost::apply_visitor(sourceNodeFetcher, b)
            );

            auto aJunctionChild = junctionInfo.firstPath.back();
            auto bJunctionChild = junctionInfo.secondPath.back();

            /* Do not use _omniDirectionalRank for root-level ranking, it should
             * only establish differences within branches
             */
            if(junctionInfo.junction != 0) {
              auto relativeRank = _omniDirectionalRank(
                junctionInfo.junction,
                TemplateMagic::setDifference(
                  _children(junctionInfo.junction),
                  std::set<TreeVertexIndex> {
                    aJunctionChild,
                    bJunctionChild
                  }
                )
              );

              /* relativeRank can only have sizes 1 or 2, where size 1 means
               * that no difference was found
               *
               * TODO this may be an accidental inversion of priority
               */
              if(relativeRank.size() == 2) {
                if(relativeRank.front().front() == aJunctionChild) {
                  relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
                } else {
                  relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
                }
              }
            }
          }
        );
      }

      // Pick the representative stereodescriptor for each branch!
      auto undecidedStereocenterSets = relativeOrders.at(branchIndex).getUndecidedSets();
      if(undecidedStereocenterSets.size() == 0) {
        // 0 - No stereocenters in the branch, so none are representative
        representativeStereodescriptors[branchIndex] = {};
      } else if(undecidedStereocenterSets.back().size() == 1) {
        // 1 - The solitary highest ranked stereogenic group
        representativeStereodescriptors[branchIndex] = {
          undecidedStereocenterSets.back().front()
        };
      } else {
        // 2 - The stereodescriptor(s) occurring more often than all others
        auto groupedByStringRep = TemplateMagic::groupByMapping(
          stereocenterMap.at(branchIndex),
          [&](const auto& variantType) -> std::string {
            return boost::apply_visitor(stringRepFetcher, variantType);
          }
        );

        auto maxSize = TemplateMagic::accumulate(
          groupedByStringRep,
          0u,
          [](const unsigned& maxSize, const auto& stringGroup) -> unsigned {
            if(stringGroup.size() > maxSize) {
              return stringGroup.size();
            }

            return maxSize;
          }
        );

        representativeStereodescriptors[branchIndex] = {};

        // All stereocenter string groups with maximum size are representative
        for(const auto& stringGroup : groupedByStringRep) {
          if(stringGroup.size() == maxSize) {
            representativeStereodescriptors.at(branchIndex).insert(
              stringGroup.begin(),
              stringGroup.end()
            );
          }
        }
      }
    }

    auto undecidedBranchSets = orderingHelper.getUndecidedSets();

    VariantLikePair variantLikeComparator {*this};

    for(const auto& undecidedSet : undecidedBranchSets) {
      TemplateMagic::forAllPairs(
        undecidedSet,
        [&](const auto& branchA, const auto& branchB) {
          // Precedence via amount of representative stereodescriptors
          if(
            representativeStereodescriptors.at(branchA).size() 
            < representativeStereodescriptors.at(branchB).size() 
          ) {
            orderingHelper.addLessThanRelationship(branchA, branchB);
          } else if(
            representativeStereodescriptors.at(branchB).size() 
            < representativeStereodescriptors.at(branchA).size() 
          ) {
            orderingHelper.addLessThanRelationship(branchB, branchA);
          } else {
            // Compare by sequential pairing
            /* What is considered a like pair?
             *
             * IUPAC says any pairs within the sets {R, M, Z} or {S, P, E}
             * are considered like-pairs, such as RR, RZ, SP, or ES. This
             * is a little more difficult for us, since any more expansive
             * stereodescriptors where choice is non-binary are not reducible
             * to two sets. Retaining the correct behavior for tetrahedral
             * stereodescriptors is the guiding principle here.
             *
             * So, we require that CN-tetrahedral-2-0 and EZ-2-0 are
             * considered like and so are CN-tetrahedral-2-1 and EZ-2-1. 
             *
             * We go about this by deciding that two stereodescriptors compare
             * like if they have the same number of assignments and have the
             * same assignment index.
             *
             * This is convenient, but means we have to ensure assignments
             * are canonical, i.e. a correspondence with R/S makes sense.
             * TODO
             */

            auto branchAOrders = relativeOrders.at(branchA).getSets();
            auto branchBOrders = relativeOrders.at(branchB).getSets();

            // Go through the ranked stereocenter groups in DESC order!
            auto branchAStereocenterGroupIter = branchAOrders.rbegin();
            auto branchBStereocenterGroupIter = branchBOrders.rbegin();

            while(
              branchAStereocenterGroupIter != branchAOrders.rend()
              && branchBStereocenterGroupIter != branchBOrders.rend()
            ) {
              unsigned ABranchLikePairs = 0, BBranchLikePairs = 0;

              // Count A-branch like pairs
              TemplateMagic::forAllPairs(
                *branchAStereocenterGroupIter,
                representativeStereodescriptors.at(branchA),
                [&](const auto& variantA, const auto& variantB) {
                  if(
                    boost::apply_visitor(
                      variantLikeComparator,
                      variantA,
                      variantB
                    )
                  ) {
                    ABranchLikePairs += 1;
                  }
                }
              );

              // Count B-branch like pairs
              TemplateMagic::forAllPairs(
                *branchBStereocenterGroupIter,
                representativeStereodescriptors.at(branchB),
                [&](const auto& variantA, const auto& variantB) {
                  if(
                    boost::apply_visitor(
                      variantLikeComparator,
                      variantA,
                      variantB
                    )
                  ) {
                    BBranchLikePairs += 1;
                  }
                }
              );

              if(ABranchLikePairs < BBranchLikePairs) {
                orderingHelper.addLessThanRelationship(branchA, branchB);
                break;
              } else if(BBranchLikePairs < ABranchLikePairs) {
                orderingHelper.addLessThanRelationship(branchB, branchA);
                break;
              }

              ++branchAStereocenterGroupIter;
              ++branchBStereocenterGroupIter;
            }
          }
        }
      );
    }

#ifndef NDEBUG
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "Sets post sequence rule 4B: {" 
      << TemplateMagic::condenseIterable(
        TemplateMagic::map(
          orderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
#endif

    // Is Sequence Rule 4, part B enough?
    if(orderingHelper.isTotallyOrdered()) {
      return fetchResult();
    }

    // TODO 4C Pseudo-asymmetries, possibly via already-gathered information?

  } // End sequence rule 4 local scope

  // Was sequence rule 4 enough?
  if(orderingHelper.isTotallyOrdered()) {
    return fetchResult();
  }

  /* Sequence rule 5: 
   * - Atom or group with {R, M, Z} precedes {S, P, E}
   */
  { // Sequence rule 5 local scope
    // Mixed BFS with specific multiset comparator, like 4A
    using VariantType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

    std::map<
      TreeVertexIndex,
      std::multiset<VariantType, SequenceRuleFiveVariantComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    // Initialize the BFS state
    for(const auto& undecidedSet : orderingHelper.getUndecidedSets()) {
      for(const auto& undecidedBranch : undecidedSet) {
        comparisonSets.emplace(
          undecidedBranch,
          *this
        );

        comparisonSets.at(undecidedBranch).emplace(undecidedBranch);
        comparisonSets.at(undecidedBranch).emplace(
          boost::edge(0, undecidedBranch, _tree).first
        );

        seeds[undecidedBranch] = {undecidedBranch};
      }
    }

    // Undecided sets not up-to-date from previous sequence rule
    undecidedSets = orderingHelper.getUndecidedSets();

    // First comparison
    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    // Recalculate undecided sets
    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 5", {0}, visitedVertices(seeds, undecidedSets)),
        _makeGraph("Sequence rule 5 multisets", 0, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a full BFS Step on all undecided set seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {
            for( // Out-edges
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdge = *outIterPair.first;

              auto edgeTarget = boost::target(outEdge, _tree);

              if(
                comparisonSets.at(undecidedTreeIndex).count(edgeTarget) == 0
                && !_tree[edgeTarget].isDuplicate
              ) {
                comparisonSets.at(undecidedTreeIndex).emplace(edgeTarget);
                comparisonSets.at(undecidedTreeIndex).emplace(outEdge);

                newSeeds.insert(edgeTarget);
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 5", {0}, visitedVertices(seeds, undecidedSets)),
          _makeGraph("Sequence rule 5 multisets", 0, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }
  } // End sequence rule 5 local scope

#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "Sets post sequence rule 5: {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

  // Exhausted sequence rules, anything undecided is now equal
  return fetchResult();
}

std::vector<
  std::vector<RankingTree::TreeVertexIndex>
> RankingTree::_omniDirectionalRank(
  const RankingTree::TreeVertexIndex& sourceIndex,
  const std::set<RankingTree::TreeVertexIndex>& excludeIndices
) const {
  /* Sequence rule 1 */
  auto sourceAdjacents = TemplateMagic::moveIf(
    _adjacents(sourceIndex),
    [&](const auto& nodeIndex) -> bool {
      // In case we explicitly exclude them, immediately discard
      if(excludeIndices.count(nodeIndex) > 0) {
        return false;
      }

      // If they are duplicate, we have some work to do
      if(_tree[nodeIndex].isDuplicate) {
        /* We need to distinguish between duplicate atoms that are the result
         * of multiple bond splits and cycle closures. Cycle closures need
         * to be kept, while multiple bond splits need to be removed. Detecting
         * the difference can be done by checking for a non-duplicate node
         * with the same molIndex adjacent to sourceIndex
         */
        auto inEdgesIterPair = boost::in_edges(sourceIndex, _tree);
        while(inEdgesIterPair.first != inEdgesIterPair.second) {
          auto adjacentIndex = boost::source(*inEdgesIterPair.first, _tree);
          if(
            _tree[adjacentIndex].molIndex == _tree[nodeIndex].molIndex
            && !_tree[adjacentIndex].isDuplicate
          ) {
            return false;
          }

          ++inEdgesIterPair.first;
        }

        auto outEdgesIterPair = boost::out_edges(sourceIndex, _tree);
        while(outEdgesIterPair.first != outEdgesIterPair.second) {
          auto adjacentIndex = boost::target(*outEdgesIterPair.first, _tree);

          if(
            _tree[adjacentIndex].molIndex == _tree[nodeIndex].molIndex
            && !_tree[adjacentIndex].isDuplicate
          ) {
            return false;
          }

          ++outEdgesIterPair.first;
        }

        /* If no non-duplicate vertex for this duplicate molIndex is immediately
         * adjacent, it is not the result of a multiple-bond split, but a cycle
         * closure, and we have to keep it
         */
        return true;
      }

      // Keep the vertex in all other cases
      return true;
    }
  );

  OrderDiscoveryHelper<TreeVertexIndex> orderingHelper {sourceAdjacents};

#ifndef NEDEBUG
Log::log(Log::Particulars::RankingTreeDebugInfo)
  << "  Preparing sequence rule 3, ranking substituents of tree index " 
  << sourceIndex 
  <<  " {" 
  << TemplateMagic::condenseIterable(
    TemplateMagic::map(
      orderingHelper.getSets(),
      [](const auto& indexSet) -> std::string {
        return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
      }
    )
  ) << "}\n";
#endif


  auto undecidedSets = orderingHelper.getUndecidedSets();

  { // local scope to avoid namespace pollution
    std::set<TreeVertexIndex> visitedVertices {sourceIndex};

    std::map<
      TreeVertexIndex,
      std::multiset<TreeVertexIndex, SequenceRuleOneVertexComparator>
    > comparisonSets;

    for(const auto& adjacent : sourceAdjacents) {
      comparisonSets.emplace(
        adjacent,
        *this
      );

      comparisonSets.at(adjacent).insert(adjacent);

      visitedVertices.insert(adjacent);
    }

    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    for(const auto& adjacentIndex : sourceAdjacents) {
      auto nextShell = _adjacents(adjacentIndex);

      // Valid seeds need to not have been encountered yet
      std::set<TreeVertexIndex> seedIndices;

      for(const auto& newSeed : nextShell) {
        if(visitedVertices.count(newSeed) == 0) {
          seedIndices.insert(newSeed);
          visitedVertices.insert(newSeed);
        }
      }

      seeds.emplace(
        adjacentIndex,
        std::move(seedIndices)
      );
    }

    // Perform the first comparison
    TemplateMagic::forAllPairs(
      sourceAdjacents,
      [&](const TreeVertexIndex& a, const TreeVertexIndex& b) {
        if(_multisetCompare(comparisonSets.at(a), comparisonSets.at(b))) {
          orderingHelper.addLessThanRelationship(a, b);
        } else if(_multisetCompare(comparisonSets.at(b), comparisonSets.at(a))) {
          orderingHelper.addLessThanRelationship(b, a);
        }
      }
    );

    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 3 prep", {0}),
        dumpGraphviz("_omni Sequence rule 1", {sourceIndex}, visitedVertices),
        _makeGraph("_omni Sequence rule 1 multisets", sourceIndex, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a BFS step on all seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {
            // Add the seed to the indices used for comparison
            comparisonSets.at(undecidedTreeIndex).insert(seed);

            // Add its un-encountered adjacents to the new seeds
            auto nextAdjacents = _adjacents(seed);

            for(const auto& newSeed : nextAdjacents) {
              if(visitedVertices.count(newSeed) == 0) {
                newSeeds.insert(newSeed);
                visitedVertices.insert(newSeed);
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 3 prep", {0}),
          dumpGraphviz("_omni Sequence rule 1", {sourceIndex}, visitedVertices),
          _makeGraph("_omni Sequence rule 1 multisets", sourceIndex, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }
  } // End sequence rule 1 scope

#ifndef NDEBUG
Log::log(Log::Particulars::RankingTreeDebugInfo)
  << "  Sets post sequence rule 1: {" 
  << TemplateMagic::condenseIterable(
    TemplateMagic::map(
      orderingHelper.getSets(),
      [](const auto& indexSet) -> std::string {
        return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
      }
    )
  ) << "}\n";
#endif

  // Is Sequence Rule 1 enough?
  if(orderingHelper.isTotallyOrdered()) {
    // No conversion of indices in _omniDirectionalRank()!
    return orderingHelper.getSets();
  }

  /* Sequence rule 2
   * - A node with higher atomic mass precedes ones with lower atomic mass
   *
   * NOTE: We skip this sequence rule for two reasons
   * - There is currently no way to represent isotopes in the library
   * - This rule may come to bear when disconnecting heteroaromatic cycles,
   *   where all aromatic atoms are given a duplicate atom with a mass equal
   *   to the average of what it would have if the double bonds were located
   *   at each of the possible positions. Although differentiating these
   *   cases would be good, there is currently no code to detect aromaticity
   *   nor permute all possible arrangements of double bonds in their
   *   subgraphs.
   */

  /* Sequence rule 3: double bond configurations
   * - Z > E > unassigned > non-stereogenic 
   *   (seqCis > seqTrans > non-stereogenic)
   */
  { // Sequence rule 3 local scope
    std::set<TreeVertexIndex> visitedVertices {sourceIndex};

    // Edge BFS outwards from each
    std::map<
      TreeVertexIndex,
      std::multiset<TreeEdgeIndex, SequenceRuleThreeEdgeComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    for(const auto& undecidedSet : undecidedSets) {
      for(const auto& undecidedVertex : undecidedSet) {
        comparisonSets.emplace(
          undecidedVertex,
          *this
        );

        auto forwardEdge = boost::edge(sourceIndex, undecidedVertex, _tree);
        auto backwardEdge = boost::edge(undecidedVertex, sourceIndex, _tree);

        comparisonSets.at(undecidedVertex).insert(
          forwardEdge.second
          ? forwardEdge.first
          : backwardEdge.first
        );

        visitedVertices.insert(undecidedVertex);

        seeds[undecidedVertex] = {undecidedVertex};
      }
    }

    // Perform first comparison
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

    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 3 prep", {0}),
        dumpGraphviz("_omni Sequence rule 3", {sourceIndex}, visitedVertices),
        _makeGraph("_omni Sequence rule 3 multisets", sourceIndex, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a BFS step on all seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {
            /* Add all unencountered adjacent edges on the seed to the
             * comparison multiset
             */

            for(
              auto inIterPair = boost::in_edges(seed, _tree);
              inIterPair.first != inIterPair.second;
              ++inIterPair.first
            ) {
              const auto& inEdgeIndex = *inIterPair.first;

              auto edgeSourceIndex = boost::source(inEdgeIndex, _tree);

              // Ensure it's not been visited before
              if(visitedVertices.count(edgeSourceIndex) == 0) {
                comparisonSets.at(undecidedTreeIndex).insert(inEdgeIndex);

                // Source nodes of an in-edge cannot be terminal, so insert
                newSeeds.insert(edgeSourceIndex);
                visitedVertices.insert(edgeSourceIndex);
              }
            }

            for(
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdgeIndex = *outIterPair.first;

              comparisonSets.at(undecidedTreeIndex).insert(outEdgeIndex);
              
              auto targetIndex = boost::target(outEdgeIndex, _tree);

              // Ensure it's not been visited before
              if(visitedVertices.count(targetIndex) == 0) {
                // Add the target node as a seed only if it's non-terminal
                if(boost::out_degree(targetIndex, _tree) > 0) {
                  newSeeds.insert(targetIndex);
                  visitedVertices.insert(targetIndex);
                }
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 3 prep", {0}),
          dumpGraphviz("_omni Sequence rule 3", {sourceIndex}, visitedVertices),
          _makeGraph("_omni Sequence rule 3 multisets", sourceIndex, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }
  } // End sequence rule 3 local scope

#ifndef NDEBUG
Log::log(Log::Particulars::RankingTreeDebugInfo)
  << "  Sets post sequence rule 3: {" 
  << TemplateMagic::condenseIterable(
    TemplateMagic::map(
      orderingHelper.getSets(),
      [](const auto& indexSet) -> std::string {
        return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
      }
    )
  ) << "}\n";
#endif

  // Is sequence rule 3 enough?
  if(orderingHelper.isTotallyOrdered()) {
    // No conversion of indices in _omniDirectionalRank()!
    return orderingHelper.getSets();
  }

  /* Sequence rule 4: Other configurations (RS, MP, EZ)
   * - stereogenic > pseudostereogenic > non-stereogenic
   * - like pairs of descriptors over unlike pairs of descriptors
   *   set 1: {R, M, Z}, set 2: {S, P, E}, like within a set, unlike crossover
   *
   *   Procedure:
   *
   *   - Create a set of all stereogenic groups within the branches to rank
   *   - Establish relative rank by application of sequence rules. 
   *
   *     Can probably use _omniDirectionalRank recursively at the first
   *     junction between competing stereocenters on the respective branch
   *     indices to find out which precedes which, provided both are at the
   *     same depth from root.
   *
   *     (WARNING: This may be tricky!)
   *
   *   - Choose a representative stereodescriptor for each branch according to
   *     the following rules:
   *
   *     - The solitarily highest ranked stereogenic group
   *     - The stereodescriptor occurring more often than all others in 
   *       the set of stereodescriptors
   *     - All most-often-occurring stereodescriptors
   *
   *   - Precedence:
   *     - The branch with fewer representative stereodescriptors has priority
   *     - Go through both lists of ranked stereocenters simultaneously,
   *       proceeding from highest rank to lowest, pairing each occurring
   *       descriptor with the reference descriptor.
   *
   *       like pairs precede unlike pairs of stereodescriptors
   *
   *
   * - (Pseudoasymmetries) r precedes s and m precedes p
   */
  { // Sequence rule 4 local scope
    /* First part of rule: A regular omni-directional BFS, inserting both
     * edges and node indices into the multiset! Perhaps avoiding the
     * insertion of uninstantiated edges and vertices entirely?
     */

    using VariantType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

    std::set<TreeVertexIndex> visitedVertices {sourceIndex};

    std::map<
      TreeVertexIndex,
      std::multiset<VariantType, SequenceRuleFourVariantComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    // Initialize the BFS state
    for(const auto& undecidedSet : orderingHelper.getUndecidedSets()) {
      for(const auto& undecidedBranch : undecidedSet) {
        comparisonSets.emplace(
          undecidedBranch,
          *this
        );

        auto forwardEdge = boost::edge(sourceIndex, undecidedBranch, _tree);
        auto backwardEdge = boost::edge(undecidedBranch, sourceIndex, _tree);

        comparisonSets.at(undecidedBranch).emplace(undecidedBranch);
        comparisonSets.at(undecidedBranch).emplace(
          forwardEdge.second
          ? forwardEdge.first
          : backwardEdge.first
        );

        seeds[undecidedBranch] = {undecidedBranch};

        visitedVertices.insert(undecidedBranch);
      }
    }

    // First comparison (undecidedSets is up-to-date from previous BFS)
    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 3 prep", {0}),
        dumpGraphviz("_omni Sequence rule 4A", {sourceIndex}, visitedVertices),
        _makeGraph("_omni Sequence rule 4A multisets", sourceIndex, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a full BFS Step on all undecided set seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {

            for( // In-edges
              auto inIterPair = boost::in_edges(seed, _tree);
              inIterPair.first != inIterPair.second;
              ++inIterPair.first
            ) {
              const auto& inEdge = *inIterPair.first;

              auto edgeSource = boost::source(inEdge, _tree);

              // Check if already placed
              if(visitedVertices.count(edgeSource) == 0) {
                comparisonSets.at(undecidedTreeIndex).emplace(edgeSource);
                comparisonSets.at(undecidedTreeIndex).emplace(inEdge);

                newSeeds.insert(edgeSource);
                visitedVertices.insert(edgeSource);
              }
            }

            for( // Out-edges
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdge = *outIterPair.first;

              auto edgeTarget = boost::target(outEdge, _tree);

              // Check if already placed and non-duplicate
              if(
                visitedVertices.count(edgeTarget) == 0
                && !_tree[edgeTarget].isDuplicate
              ) {
                comparisonSets.at(undecidedTreeIndex).emplace(edgeTarget);
                comparisonSets.at(undecidedTreeIndex).emplace(outEdge);

                newSeeds.insert(edgeTarget);
                visitedVertices.insert(edgeTarget);
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 3 prep", {0}),
          dumpGraphviz("_omni Sequence rule 4A", {sourceIndex}, visitedVertices),
          _makeGraph("_omni Sequence rule 4A multisets", sourceIndex, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif
    }

#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "  Sets post sequence rule 4A: {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

    // Is Sequence Rule 4, part A enough?
    if(orderingHelper.isTotallyOrdered()) {
      // No conversion of indices in _omniDirectionalRank()!
      return orderingHelper.getSets();
    }

    /* Second part: Choosing representational stereodescriptors for all
     * branches and pairing in sequence of priority for branches.
     */
    std::map<
      TreeVertexIndex,
      std::set<VariantType>
    > stereocenterMap;

    VariantHasInstantiatedStereocenter isInstantiatedChecker {*this};

    // Copy all those variants from part A that are actually instantiated
    for(const auto& undecidedSet : orderingHelper.getUndecidedSets()) {
      for(const auto& undecidedBranch : undecidedSet) {
        stereocenterMap[undecidedBranch] = {};

        for(const auto& variantValue : comparisonSets.at(undecidedBranch)) {
          if(boost::apply_visitor(isInstantiatedChecker, variantValue)) {
            stereocenterMap.at(undecidedBranch).insert(variantValue);
          }
        }
      }
    }

    std::map<
      TreeVertexIndex,
      OrderDiscoveryHelper<VariantType>
    > relativeOrders;

    VariantDepth depthFetcher {*this};
    VariantSourceNode sourceNodeFetcher {*this};
    VariantStereocenterStringRepresentation stringRepFetcher {*this};

    std::map<
      TreeVertexIndex,
      std::set<VariantType>
    > representativeStereodescriptors;

    for(const auto& mapIterPair : stereocenterMap) {
      const auto& branchIndex = mapIterPair.first;
      const auto& variantSet = mapIterPair.second;

      // Instantiate the order discovery helpers
      relativeOrders.emplace(
        branchIndex,
        variantSet
      );

      // Compare based on depth
      TemplateMagic::forAllPairs(
        variantSet,
        [&](const auto& a, const auto& b) {
          auto aDepth = boost::apply_visitor(depthFetcher, a);
          auto bDepth = boost::apply_visitor(depthFetcher, b);

          if(aDepth < bDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
          } else if(bDepth < aDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
          }
        }
      );

      /* Try to resolve undecided sets via ranking downwards at junction of
       * pairs
       */
      for(
        const auto& undecidedSet 
        : relativeOrders.at(branchIndex).getUndecidedSets()
      ) {
        TemplateMagic::forAllPairs(
          undecidedSet,
          [&](const auto& a, const auto& b) {
            auto junctionInfo = _junction(
              boost::apply_visitor(sourceNodeFetcher, a),
              boost::apply_visitor(sourceNodeFetcher, b)
            );

            auto aJunctionChild = junctionInfo.firstPath.back();
            auto bJunctionChild = junctionInfo.secondPath.back();

            /* Do not use _omniDirectionalRank for root-level ranking, it should
             * only establish differences within branches here
             */
            if(junctionInfo.junction != 0) {
              auto relativeRank = _omniDirectionalRank(
                junctionInfo.junction,
                TemplateMagic::setDifference(
                  _children(junctionInfo.junction),
                  std::set<TreeVertexIndex> {
                    aJunctionChild,
                    bJunctionChild
                  }
                )
              );

              /* relativeRank can only have sizes 1 or 2, 1 meaning that no
               * difference was found
               */
              if(relativeRank.size() == 2) {
                if(relativeRank.front().front() == aJunctionChild) {
                  relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
                } else {
                  relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
                }
              }
            }
          }
        );
      }

      // Pick the representative stereodescriptor for each branch!
      auto undecidedStereocenterSets = relativeOrders.at(branchIndex).getUndecidedSets();

      if(undecidedStereocenterSets.size() == 0) {
        // 0 - No stereocenters in the branch, so none are representative
        representativeStereodescriptors[branchIndex] = {};
      } else if(undecidedStereocenterSets.back().size() == 1) {
        // 1 - The solitary highest ranked stereogenic group
        representativeStereodescriptors[branchIndex] = {
          undecidedStereocenterSets.back().front()
        };
      } else {
        // 2 - The stereodescriptor occurring more often than all others
        auto groupedByStringRep = TemplateMagic::groupByMapping(
          stereocenterMap.at(branchIndex),
          [&](const auto& variantType) -> std::string {
            return boost::apply_visitor(stringRepFetcher, variantType);
          }
        );

        auto maxSize = TemplateMagic::accumulate(
          groupedByStringRep,
          0u,
          [](const unsigned& maxSize, const auto& stringGroup) -> unsigned {
            if(stringGroup.size() > maxSize) {
              return stringGroup.size();
            }

            return maxSize;
          }
        );

        representativeStereodescriptors[branchIndex] = {};

        // All stereocenter string groups with maximum size are representative
        for(const auto& stringGroup : groupedByStringRep) {
          if(stringGroup.size() == maxSize) {
            representativeStereodescriptors.at(branchIndex).insert(
              stringGroup.begin(),
              stringGroup.end()
            );
          }
        }
      }
    }

    auto undecidedBranchSets = orderingHelper.getUndecidedSets();

    VariantLikePair variantLikeComparator {*this};

    for(const auto& undecidedSet : undecidedBranchSets) {
      TemplateMagic::forAllPairs(
        undecidedSet,
        [&](const auto& branchA, const auto& branchB) {
          // Precedence via amount of representative stereodescriptors
          if(
            representativeStereodescriptors.at(branchA).size() 
            < representativeStereodescriptors.at(branchB).size() 
          ) {
            orderingHelper.addLessThanRelationship(branchA, branchB);
          } else if(
            representativeStereodescriptors.at(branchB).size() 
            < representativeStereodescriptors.at(branchA).size() 
          ) {
            orderingHelper.addLessThanRelationship(branchB, branchA);
          } else {
            // Compare by sequential pairing
            /* What is considered a like pair?
             *
             * IUPAC says any pairs within the sets {R, M, Z} or {S, P, E}
             * are considered like-pairs, such as RR, RZ, SP, or ES. This
             * is a little more difficult for us, since any more expansive
             * stereodescriptors where choice is non-binary are not reducible
             * to two sets. Retaining the correct behavior for tetrahedral
             * stereodescriptors is the guiding principle here.
             *
             * So, we require that CN-tetrahedral-2-0 and EZ-2-0 are
             * considered like and so are CN-tetrahedral-2-1 and EZ-2-1. 
             *
             * We go about this by deciding that two stereodescriptors compare
             * like if they have the same number of assignments and have the
             * same assignment index.
             *
             * This is convenient, but means we have to ensure assignments
             * are canonical, i.e. a correspondence with R/S makes sense.
             * TODO
             */

            auto branchAOrders = relativeOrders.at(branchA).getSets();
            auto branchBOrders = relativeOrders.at(branchB).getSets();

            // Go through the ranked stereocenter groups in DESC order!
            auto branchAStereocenterGroupIter = branchAOrders.rbegin();
            auto branchBStereocenterGroupIter = branchBOrders.rbegin();

            while(
              branchAStereocenterGroupIter != branchAOrders.rend()
              && branchBStereocenterGroupIter != branchBOrders.rend()
            ) {
              unsigned ABranchLikePairs = 0, BBranchLikePairs = 0;

              // Count A-branch like pairs
              TemplateMagic::forAllPairs(
                *branchAStereocenterGroupIter,
                representativeStereodescriptors.at(branchA),
                [&](const auto& variantA, const auto& variantB) {
                  if(
                    boost::apply_visitor(
                      variantLikeComparator,
                      variantA,
                      variantB
                    )
                  ) {
                    ABranchLikePairs += 1;
                  }
                }
              );

              // Count B-branch like pairs
              TemplateMagic::forAllPairs(
                *branchBStereocenterGroupIter,
                representativeStereodescriptors.at(branchB),
                [&](const auto& variantA, const auto& variantB) {
                  if(
                    boost::apply_visitor(
                      variantLikeComparator,
                      variantA,
                      variantB
                    )
                  ) {
                    BBranchLikePairs += 1;
                  }
                }
              );

              if(ABranchLikePairs < BBranchLikePairs) {
                orderingHelper.addLessThanRelationship(branchA, branchB);
                break;
              } else if(BBranchLikePairs < ABranchLikePairs) {
                orderingHelper.addLessThanRelationship(branchB, branchA);
                break;
              }

              ++branchAStereocenterGroupIter;
              ++branchBStereocenterGroupIter;
            }
          }
        }
      );
    }

#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "  Sets post sequence rule 4B: {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

    // Is Sequence Rule 4, part B enough?
    if(orderingHelper.isTotallyOrdered()) {
      // No conversion of indices in _omniDirectionalRank()!
      return orderingHelper.getSets();
    }

    // TODO 4C Pseudo-asymmetries, possibly via already-gathered information?

  } // End sequence rule 4 local scope

  /* Sequence rule 5: 
   * - Atom or group with {R, M, Z} precedes {S, P, E}
   */
  { // Sequence rule 5 local scope
    // Mixed BFS with specific multiset comparator, like 4A
    using VariantType = boost::variant<TreeVertexIndex, TreeEdgeIndex>;

    std::set<TreeVertexIndex> visitedVertices {sourceIndex};

    std::map<
      TreeVertexIndex,
      std::multiset<VariantType, SequenceRuleFiveVariantComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::set<TreeVertexIndex>
    > seeds;

    // Initialize the BFS state
    for(const auto& undecidedSet : orderingHelper.getUndecidedSets()) {
      for(const auto& undecidedBranch : undecidedSet) {
        comparisonSets.emplace(
          undecidedBranch,
          *this
        );

        auto forwardEdge = boost::edge(sourceIndex, undecidedBranch, _tree);
        auto backwardEdge = boost::edge(undecidedBranch, sourceIndex, _tree);
        auto edgeIndex = forwardEdge.second
          ? forwardEdge.first
          : backwardEdge.first;

        comparisonSets.at(undecidedBranch).emplace(undecidedBranch);
        comparisonSets.at(undecidedBranch).emplace(edgeIndex);

        seeds[undecidedBranch] = {undecidedBranch};

        visitedVertices.insert(undecidedBranch);
      }
    }

    undecidedSets = orderingHelper.getUndecidedSets();

    // First comparison
    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    // Recalculate undecided sets
    undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz("Sequence rule 3 prep", {0}),
        dumpGraphviz("_omni Sequence rule 5", {sourceIndex}, visitedVertices),
        _makeGraph("_omni Sequence rule 5 multisets", sourceIndex, comparisonSets, undecidedSets),
        orderingHelper.dumpGraphviz()
      });
    }
#endif

    while(undecidedSets.size() > 0 && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a full BFS Step on all undecided set seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::set<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {

            for( // In-edges
              auto inIterPair = boost::in_edges(seed, _tree);
              inIterPair.first != inIterPair.second;
              ++inIterPair.first
            ) {
              const auto& inEdge = *inIterPair.first;

              auto edgeSource = boost::source(inEdge, _tree);

              if(visitedVertices.count(edgeSource) == 0) {

                comparisonSets.at(undecidedTreeIndex).emplace(edgeSource);
                comparisonSets.at(undecidedTreeIndex).emplace(inEdge);

                newSeeds.insert(edgeSource);
                visitedVertices.insert(edgeSource);
              }
            }

            for( // Out-edges
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdge = *outIterPair.first;

              auto edgeTarget = boost::target(outEdge, _tree);

              if(
                visitedVertices.count(edgeTarget) == 0
                && !_tree[edgeTarget].isDuplicate
              ) {
                comparisonSets.at(undecidedTreeIndex).emplace(edgeTarget);
                comparisonSets.at(undecidedTreeIndex).emplace(outEdge);

                newSeeds.insert(edgeTarget);
                visitedVertices.insert(edgeTarget);
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

      // Recalculate the undecided sets
      undecidedSets = orderingHelper.getUndecidedSets();

#ifndef NDEBUG
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptMolGraph(_moleculeRef.dumpGraphviz()),
          dumpGraphviz("Sequence rule 3 prep", {0}),
          dumpGraphviz("_omni Sequence rule 5", {sourceIndex}, visitedVertices),
          _makeGraph("_omni Sequence rule 5 multisets", sourceIndex, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
#endif

    }
  } // End sequence rule 5 local scope

#ifndef NDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "  Sets post sequence rule 5: {" 
    << TemplateMagic::condenseIterable(
      TemplateMagic::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + TemplateMagic::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "}\n";
#endif

  // Exhausted sequence rules, return the sets
  return orderingHelper.getSets();
}

void RankingTree::_DFSExpand(
  const RankingTree::TreeVertexIndex& index,
  const std::set<AtomIndexType>& molIndicesInBranch
) {
  std::set<AtomIndexType> treeOutAdjacencies;
  TreeGraphType::out_edge_iterator iter, end;
  std::tie(iter, end) = boost::out_edges(index, _tree);
  while(iter != end) {
    treeOutAdjacencies.insert(
      _tree[
        boost::target(*iter, _tree)
      ].molIndex
    );

    ++iter;
  }

  for(
    const auto& molAdjacentIndex 
    : _moleculeRef.iterateAdjacencies(_tree[index].molIndex)
  ) {
    if(treeOutAdjacencies.count(molAdjacentIndex) != 0) {
      continue;
    }

    if(molIndicesInBranch.count(molAdjacentIndex) == 0) {
      auto newIndex = boost::add_vertex(_tree);
      _tree[newIndex].molIndex = molAdjacentIndex;
      _tree[newIndex].isDuplicate = false;

      boost::add_edge(index, newIndex, _tree);
      
      // Need to add duplicates!
      _addBondOrderDuplicates(index, newIndex);

      // Copy the indices and insert the newest addition
      auto molIndicesInBranchCopy = molIndicesInBranch;
      molIndicesInBranchCopy.insert(molAdjacentIndex);

      _DFSExpand(newIndex, molIndicesInBranchCopy);
    } else if(molAdjacentIndex != _tree[_parent(index)].molIndex)  {
      auto newIndex = boost::add_vertex(_tree);
      _tree[newIndex].molIndex = molAdjacentIndex;
      _tree[newIndex].isDuplicate = true;

      boost::add_edge(index, newIndex, _tree);

      // Terminate DFS
    }
  }
}

template<>
std::string RankingTree::toString(const TreeVertexIndex& vertexIndex) const {
  return std::to_string(vertexIndex);
}

template<>
std::string RankingTree::toString(const TreeEdgeIndex& edgeIndex) const {
  return (
    std::to_string(
      boost::source(edgeIndex, _tree)
    ) + "->"s + std::to_string(
      boost::target(edgeIndex, _tree)
    )
  );
}

template<>
std::string RankingTree::toString(const boost::variant<TreeVertexIndex, TreeEdgeIndex>& variant) const {
  if(variant.which() == 0) {
    return toString<TreeVertexIndex>(boost::get<TreeVertexIndex>(variant));
  } 
  
  return toString<TreeEdgeIndex>(boost::get<TreeEdgeIndex>(variant));
}

//! Returns the parent of a node. Fails if called on the root!
RankingTree::TreeVertexIndex RankingTree::_parent(const RankingTree::TreeVertexIndex& index) const {
  assert(index != 0);

  // All nodes in the graph must have in_degree of 1
  assert(boost::in_degree(index, _tree) == 1);

  // Follow singular in_node to source, overwrite index
  auto iterPair = boost::in_edges(index, _tree);
  return boost::source(*iterPair.first, _tree);
}

//! Returns the direct descendants of a tree node
std::set<RankingTree::TreeVertexIndex> RankingTree::_children(const RankingTree::TreeVertexIndex& index) const {
  std::set<TreeVertexIndex> children;

  TreeGraphType::out_edge_iterator iter, end;
  std::tie(iter, end) = boost::out_edges(index, _tree);

  while(iter != end) {
    children.insert(
      boost::target(*iter, _tree)
    );

    ++iter;
  }

  return children;
}

std::set<RankingTree::TreeVertexIndex> RankingTree::_adjacents(const RankingTree::TreeVertexIndex& index) const {
  std::set<TreeVertexIndex> adjacents;

  auto inIterPair = boost::in_edges(index, _tree);
  while(inIterPair.first != inIterPair.second) {
    adjacents.emplace(
      boost::source(*inIterPair.first, _tree)
    );

    ++inIterPair.first;
  }

  auto outIterPair = boost::out_edges(index, _tree);
  while(outIterPair.first != outIterPair.second) {
    adjacents.emplace(
      boost::target(*outIterPair.first, _tree)
    );

    ++outIterPair.first;
  }

  return adjacents;
}

std::set<RankingTree::TreeEdgeIndex> RankingTree::_adjacentEdges(const RankingTree::TreeVertexIndex& index) const {
  std::set<TreeEdgeIndex> edges;

  auto inIterPair = boost::in_edges(index, _tree);
  while(inIterPair.first != inIterPair.second) {
    edges.emplace(*inIterPair.first);

    ++inIterPair.first;
  }

  auto outIterPair = boost::out_edges(index, _tree);
  while(outIterPair.first != outIterPair.second) {
    edges.emplace(*outIterPair.first);

    ++outIterPair.first;
  }

  return edges;
}

unsigned RankingTree::_nonDuplicateDegree(const RankingTree::TreeVertexIndex& index) const {
  auto adjacents = _adjacents(index);

  auto numDuplicate = TemplateMagic::accumulate(
    adjacents,
    0u,
    [&](const unsigned& count, const auto& treeIndex) -> unsigned {
      return count + static_cast<unsigned>(
        _tree[treeIndex].isDuplicate
      );
    }
  );

  return adjacents.size() - numDuplicate;
}

bool RankingTree::_relevantSeeds(
  const std::map<
    RankingTree::TreeVertexIndex, 
    std::set<RankingTree::TreeVertexIndex>
  >& seeds,
  const std::vector<
    std::vector<RankingTree::TreeVertexIndex>
  >& undecidedSets
) {
  for(const auto& undecidedSet : undecidedSets) {
    for(const auto& undecidedTreeIndex : undecidedSet) {
      if(seeds.at(undecidedTreeIndex).size() > 0) {
        return true;
      }
    }
  }

  return false;
}

void RankingTree::_addBondOrderDuplicates(
  const TreeVertexIndex& treeSource,
  const TreeVertexIndex& treeTarget
) {
  const auto& molGraph = _moleculeRef.getGraph();

  auto molGraphEdge = boost::edge(
    _tree[treeSource].molIndex,
    _tree[treeTarget].molIndex,
    molGraph
  );

  /* In case the bond order is non-fractional (aromatic / eta)
   * and > 1, add duplicate atoms
   */
  assert(molGraphEdge.second);

  auto bondType = molGraph[molGraphEdge.first].bondType;

  /* If the bond is double, we must remember it for sequence rule 2.
   * Remembering these edges is way more convenient than looking for subgraphs
   * which have the trace pattern of a split multiple-bond of BO 2.
   */
  if(bondType == BondType::Double) {
    auto edgePair = boost::edge(treeSource, treeTarget, _tree);
    _doubleBondEdges.emplace(edgePair.first);
  }

  double bondOrder = Bond::bondOrderMap.at(
    static_cast<unsigned>(bondType)
  );
  unsigned integralBondOrder = static_cast<unsigned>(bondOrder);

  // Check if bond order is integral
  if(static_cast<double>(integralBondOrder) == bondOrder) {
    // Add duplicate atom to source and target integralBondOrder - 1 times
    for(unsigned N = 1; N < integralBondOrder; ++N) {
      // Add a node with molIndex molSource to treeTarget
      auto a = boost::add_vertex(_tree);
      _tree[a].molIndex = _tree[treeSource].molIndex;
      _tree[a].isDuplicate = true;
      boost::add_edge(treeTarget, a, _tree);

      // Add a node with molIndex molTarget to treeSource
      auto b = boost::add_vertex(_tree);
      _tree[b].molIndex = _tree[treeTarget].molIndex;
      _tree[b].isDuplicate = true;

      boost::add_edge(treeSource, b, _tree);
    }
  }
}

std::set<AtomIndexType> RankingTree::_molIndicesInBranch(TreeVertexIndex index) const {
  std::set<AtomIndexType> indices {_tree[index].molIndex};

  while(index != 0) {
    // All nodes in the graph must have in_degree of 1
    assert(boost::in_degree(index, _tree) == 1);

    // Follow singular in_node to source, overwrite index
    auto iterPair = boost::in_edges(index, _tree);
    index = boost::source(*iterPair.first, _tree);

    indices.insert(_tree[index].molIndex);
  }

  return indices;
}

unsigned RankingTree::_duplicateDepth(TreeVertexIndex index) const {
  // Call this only on duplicate vertices!
  assert(_tree[index].isDuplicate);

  auto duplicateMolIndex = _tree[index].molIndex;

  /* Summary:
   * - If the original atom is not found, return actual depth
   * - If the original atom is found, return that vertex's depth
   */

  unsigned depth = 0;

  while(index != 0) {
    // All nodes in the graph must have in_degree of 1
    assert(boost::in_degree(index, _tree) == 1);

    // Follow singular in_node to source, overwrite index
    auto iterPair = boost::in_edges(index, _tree);
    index = boost::source(*iterPair.first, _tree);

    ++depth;

    if(_tree[index].molIndex == duplicateMolIndex) {
      /* Reset depth. Since the loop continues until root, this returns the
       * depth of the non-duplicate node
       */
      depth = 0;
    }
  }

  return depth;
}

//! Returns the depth of a node in the tree
unsigned RankingTree::_depthOfNode(TreeVertexIndex index) const {
  /* As long as VertexListS is vecS, the first node (root) will always have
   * index zero
   */
  unsigned depth = 0;

  while(index != 0) {
    // All nodes in the graph must have in_degree of 1
    assert(boost::in_degree(index, _tree) == 1);

    // Follow singular in_node to source, overwrite index
    auto iterPair = boost::in_edges(index, _tree);
    index = boost::source(*iterPair.first, _tree);

    ++depth;
  }

  return depth;
}

//!  Returns a mixed depth measure for ranking both vertices and edges
unsigned RankingTree::_mixedDepth(const TreeVertexIndex& vertexIndex) const {
  return 2 * _depthOfNode(vertexIndex);
}

//!  Returns a mixed depth measure for ranking both vertices and edges
unsigned RankingTree::_mixedDepth(const TreeEdgeIndex& edgeIndex) const {
  return 2 * _depthOfNode(
    boost::source(edgeIndex, _tree)
  ) + 1;
}

struct RankingTree::JunctionInfo {
  TreeVertexIndex junction;

  std::vector<TreeVertexIndex> firstPath, secondPath;
};

/*!
 * Returns the deepest vertex in the tree in whose child branches both a and
 * b are located
 */
typename RankingTree::JunctionInfo RankingTree::_junction(
  const TreeVertexIndex& a,
  const TreeVertexIndex& b
) const {
  JunctionInfo data;


  /* Determine the junction vertex */
  std::set<TreeVertexIndex> aBranchIndices = _molIndicesInBranch(a);

  // By default, the junction is root
  data.junction = 0;

  // In case b is included in the path, that is the junction
  if(aBranchIndices.count(b)) {
    data.junction = b;
  } else {
    // Backtrack with b, checking at each iteration
    auto bCurrent = b;
    while(bCurrent != 0) {
      bCurrent = _parent(bCurrent);

      if(aBranchIndices.count(bCurrent)) {
        data.junction = bCurrent;
      }
    }
  }

  /* Reconstruct the paths to the junction */

  for(
    TreeVertexIndex aCurrent = a;
    a != data.junction;
    aCurrent = _parent(aCurrent)
  ) {
    data.firstPath.push_back(aCurrent);
  }

  for(
    TreeVertexIndex bCurrent = b;
    b != data.junction;
    bCurrent = _parent(bCurrent)
  ) {
    data.secondPath.push_back(bCurrent);
  }

  return data;
}

//! Returns whether a molecular graph index exists in a specific branch
bool RankingTree::_molIndexExistsInBranch(
  const AtomIndexType& molIndex,
  TreeVertexIndex treeIndex
) const {
  // Check whether this treeIndex already has the molIndex
  if(_tree[treeIndex].molIndex == molIndex) {
    return true;
  }

  while(treeIndex != 0) {
    // All nodes in the graph must have in_degree of 1
    assert(boost::in_degree(treeIndex, _tree) == 1);

    /* Move up one node in the graph. Follow singular in_node to source,
     * overwrite index
     */
    auto iterPair = boost::in_edges(treeIndex, _tree);
    treeIndex = boost::source(*iterPair.first, _tree);

    // Test this position
    if(_tree[treeIndex].molIndex == molIndex) {
      return true;
    }
  }

  return false;
}

std::string RankingTree::dumpGraphviz(
  const std::string& title,
  const std::set<TreeVertexIndex>& squareVertices,
  const std::set<TreeVertexIndex>& colorVertices,
  const std::set<TreeEdgeIndex>& colorEdges
) const {
  GraphvizWriter propertyWriter {
    *this,
    title,
    squareVertices,
    colorVertices,
    colorEdges
  };

  std::stringstream ss;

  boost::write_graphviz(
    ss,
    _tree,
    propertyWriter,
    propertyWriter,
    propertyWriter
  );

  return ss.str();
}

const typename RankingTree::TreeGraphType& RankingTree::getGraph() const {
  return _tree;
}

#ifndef NDEBUG
// Initialize the debug counter
unsigned RankingTree::_debugMessageCounter = 0;

void RankingTree::_writeGraphvizFiles(
  const std::vector<std::string>& graphvizStrings
) {
  using namespace std::string_literals;

  for(unsigned i = 0; i < graphvizStrings.size(); ++i) {
    std::string filename = (
      "ranking-tree-"s
      + std::to_string(_debugMessageCounter)
      + "-"s
      + std::to_string(i)
      + ".dot"s
    );


    std::ofstream outFile(filename);
    outFile << graphvizStrings.at(i);
    outFile.close();
  }

  ++_debugMessageCounter;
}

std::string RankingTree::_adaptMolGraph(std::string molGraph) {
  /* We have to make it a digraph for compatibility with all other graphs, so:
   * 1: "graph G {" -> "digraph G {"
   */
  molGraph.insert(0, "di");

  // 2: all "--" edges must become -> edges
  auto firstEdgePos = molGraph.find_first_of("--");
  boost::replace_all(molGraph, "--", "->");

  // 3: all edges need a dir="none" specification to avoid directed edges
  auto bracketPos = molGraph.find_first_of(']', firstEdgePos);
  while(bracketPos != std::string::npos) {
    molGraph.insert(bracketPos, R"(, dir="none")");
    bracketPos += 13;
    bracketPos = molGraph.find_first_of(']', bracketPos);
  }

  return molGraph;
}
#endif


} // namespace MoleculeManip
