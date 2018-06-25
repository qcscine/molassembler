#include "RankingTree.h"

#include "boost/algorithm/string/replace.hpp"
#include "boost/graph/breadth_first_search.hpp"
#include "Delib/ElementInfo.h"

#include "chemical_symmetries/Properties.h"
#include "temple/constexpr/ConsecutiveCompare.h"

#include "detail/MolGraphWriter.h"
#include "GraphHelpers.h"
#include "LocalGeometryModel.h"
#include "Options.h"

#include <fstream>
#include <iostream>

namespace molassembler {

// Must declare constexpr static member without definition!
constexpr decltype(RankingTree::rootIndex) RankingTree::rootIndex;

//! Helper class to write a graphviz representation of the generated tree
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
    std::string title = "",
    std::set<TreeVertexIndex> squareVertices = {},
    std::set<TreeVertexIndex> colorVertices = {},
    std::set<TreeEdgeIndex> colorEdges = {}
  ) : _baseRef(baseTree),
      _title(std::move(title)),
      _squareVertices(std::move(squareVertices)),
      _colorVertices(std::move(colorVertices)),
      _colorEdges(std::move(colorEdges))
  {}

  void operator() (std::ostream& os) const {
    os << R"(  graph [fontname = "Arial"];)" << "\n"
      << R"(  node [fontname = "Arial", shape = "circle", style = "filled"];)" << "\n"
      << R"(  edge [fontname = "Arial"];)" << "\n"
      << R"(  labelloc="t"; label=")" << _title << R"(")" << ";\n";
  }

  void operator() (std::ostream& os, const TreeVertexIndex& vertexIndex) const {
    auto symbolString = Delib::ElementInfo::symbol(
      graph::elementType(
        _baseRef._tree[vertexIndex].molIndex,
        _baseRef._graphRef
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
/*! IUPAC Sequence rule one tree vertex comparator
 *
 * Comparator class for comparing individual tree vertex indices according to
 * IUPAC Sequence rule one. Handles all combinations of duplicate and
 * non-duplicate nodes.
 *
 * @note This is a non-default-instantiable Comparator, meaning care must be
 * taken in the instantiation of the STL Container using this to avoid move
 * and copy assignment operators. You must use in-place-construction!
 */
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
        graph::elementType(
          _base._tree[a].molIndex,
          _base._graphRef
        )
      ) < static_cast<unsigned>(
        graph::elementType(
          _base._tree[b].molIndex,
          _base._graphRef
        )
      );
    }

    if(aDuplicate && bDuplicate) {
      // That vertex is smaller which is further from root (invert comparison)
      return _base._duplicateDepth(a) > _base._duplicateDepth(b);
    }

    // The case below is not needed, it falls together with equality
    /*if(!aDuplicate && bDuplicate) {
      return false;
    }*/

    return aDuplicate && !bDuplicate;
  }
};

/*! IUPAC Sequence rule three tree edge comparator
 *
 * Comparator class for comparing individual tree edges according to IUPAC
 * sequence rule three.
 *
 * @note This is a non-default-instantiable Comparator, meaning care must be
 * taken in the instantiation of the STL Container using this to avoid move
 * and copy assignment operators. You must use in-place-construction!
 */
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
     *   - Stereopermutation value (Z > E)
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
    return temple::consecutiveCompareSmaller(
      EZStereocenterB.numStereopermutations(),
      EZStereocenterA.numStereopermutations(),
      /* Mixed optional comparison (includes comparison of assignment value if
       * assigned)
       */
      EZStereocenterB.indexOfPermutation(),
      EZStereocenterA.indexOfPermutation()
    );
  }
};

/*! IUPAC Sequence rule four tree vertex and edge mixed comparator
 *
 * Comparator class for comparing both tree vertices and edges according to
 * IUPAC sequence rule four.
 *
 * @note This is a non-default-instantiable Comparator, meaning care must be
 * taken in the instantiation of the STL Container using this to avoid move
 * and copy assignment operators. You must use in-place-construction!
 */
struct RankingTree::SequenceRuleFourVariantComparator {
public:
  class VariantComparisonVisitor : boost::static_visitor<bool> {
  private:
    const RankingTree& _base;

  public:
    explicit VariantComparisonVisitor(
      const SequenceRuleFourVariantComparator& comparatorBase
    ) : _base(comparatorBase._base) {}

    /* Surprisingly, the code for homogeneous and heterogeneous comparisons is
     * completely identical, so we can abstract over the types.
     */
    template<typename T, typename U>
    bool operator() (const T& a, const U& b) const {
      auto aDepth = _base._mixedDepth(a);
      auto bDepth = _base._mixedDepth(b);

      if(aDepth < bDepth) {
        return true;
      }

      if(aDepth > bDepth) {
        return false;
      }

      const auto& aOption = _base._tree[a].stereocenterOption;
      const auto& bOption = _base._tree[b].stereocenterOption;

      // Uninstantiated stereocenters always compare false
      if(!aOption && !bOption) {
        return false;
      }

      // Instantiated precedes uninstantiated
      if(aOption && !bOption) {
        return true;
      }

      if(!aOption && bOption) {
        return false;
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
       *   false > true == 0 > 1 == false
       *   (So A < B is false, meaning stereogenic < non-stereogenic, leading
       *   to the desired ordering)
       *
       * This is valid for both CN and EZ types of stereocenters
       */
      return (
        (StereocenterA.numStereopermutations() > 1)
        > (StereocenterB.numStereopermutations() > 1)
      );
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

  bool operator () (const VariantType& a, const VariantType& b) const {
    return boost::apply_visitor(_variantComparator, a, b);
  }
};

/*! IUPAC Sequence rule five tree vertex and edge mixed comparator
 *
 * Comparator class for comparing both tree vertices and edges according to
 * IUPAC sequence rule five.
 *
 * @note This is a non-default-instantiable Comparator, meaning care must be
 * taken in the instantiation of the STL Container using this to avoid move
 * and copy assignment operators. You must use in-place-construction!
 */
struct RankingTree::SequenceRuleFiveVariantComparator {
public:
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
      if(StereocenterA.numStereopermutations() < StereocenterB.numStereopermutations()) {
        return true;
      }

      if(StereocenterB.numStereopermutations() < StereocenterA.numStereopermutations()) {
        return false;
      }

      /* In our context, sequence rule 5 directly compares the assignment
       * of the assigned stereocenters.
       */
      return StereocenterA.indexOfPermutation() < StereocenterB.indexOfPermutation();
    }

    // For different types
    template<typename T, typename U>
    bool operator() (const T& /* t */, const U& /* u */) const {
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

  bool operator () (const VariantType& a, const VariantType& b) const {
    return boost::apply_visitor(_variantComparator, a, b);
  }
};

/* Variant visitor helper classes */
//! Predicate of whether tree vertex or edge have an instantiated stereocenter
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

//! Functor calculating mixed depth (see _mixedDepth functions) of tree vertex or edge
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

//! Functor fetching the source vertex of a tree edge or vertex (identity)
class RankingTree::VariantSourceNode : boost::static_visitor<TreeVertexIndex> {
private:
  const RankingTree& _baseRef;

public:
  explicit VariantSourceNode(const RankingTree& base) : _baseRef(base) {}

  //! Returns the vertex itself
  TreeVertexIndex operator() (const TreeVertexIndex& a) const {
    return a;
  }

  //! Returns the source of the edge (vertex closer to root)
  TreeVertexIndex operator() (const TreeEdgeIndex& a) const {
    return boost::source(a, _baseRef._tree);
  }
};

//! Functor returning a string representation of a vertex or edge stereocenter
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

//! Predicate deciding IUPAC like pairing of combination of tree vertex and edge stereocenters
class RankingTree::VariantLikePair : boost::static_visitor<bool> {
private:
  const RankingTree& _baseRef;

public:
  explicit VariantLikePair(const RankingTree& base) : _baseRef(base) {}

  template<typename T, typename U>
  bool operator() (const T& a, const U& b) const {
    const auto& aOption = _baseRef._tree[a].stereocenterOption;
    const auto& bOption = _baseRef._tree[b].stereocenterOption;

    return (
      aOption
      && bOption
      && aOption.value().numStereopermutations() == bOption.value().numStereopermutations()
      && aOption.value().indexOfPermutation() == bOption.value().indexOfPermutation()
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
void RankingTree::_applySequenceRules(
  const boost::optional<AngstromWrapper>& positionsOption
) {
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

  /* A much-needed variable which was often re-declared within the local
   * scopes of every sequence rule. This allows reuse.
   */
  auto undecidedSets = _branchOrderingHelper.getUndecidedSets();

  { // Populate byDepth from active branches only
    std::set<TreeEdgeIndex> rootEdges;

    for(const auto& undecidedSet : undecidedSets) {
      for(const auto& undecidedBranch : undecidedSet) {
        rootEdges.insert(
          boost::edge(rootIndex, undecidedBranch, _tree).first
        );
      }
    }

    byDepth.emplace_back(
      std::move(rootEdges)
    );
  }

  std::set<TreeEdgeIndex> nextChildren;
  bool moreEdges;

  do {
    nextChildren.clear();

    for(const auto& edge : byDepth.back()) {
      auto outIterPair = boost::out_edges(
        boost::target(edge, _tree), _tree
      );

      while(outIterPair.first != outIterPair.second) {
        if( // Only add edges that have non-terminal targets
          boost::out_degree(
            boost::target(
              *outIterPair.first,
              _tree
            ),
            _tree
          ) > 0
        ) {
          nextChildren.insert(*outIterPair.first);
        }

        ++outIterPair.first;
      }
    }

    moreEdges = !nextChildren.empty();
    if(moreEdges) {
      byDepth.emplace_back(
        std::move(nextChildren)
      );
    }
  } while(moreEdges);

  // Remember if you found something to instantiate or not
  bool foundEZStereocenters = false;
  bool foundCNStereocenters = false;


  auto auxiliaryLigands = [](
    const RankingInformation::RankedType& rankedAtoms
  ) -> RankingInformation::LigandsType {
    RankingInformation::LigandsType ligands;

    for(const auto& equalAtomSet : rankedAtoms) {
      for(const AtomIndexType individualAtom : equalAtomSet) {
        ligands.push_back(
          std::vector<AtomIndexType> {individualAtom}
        );
      }
    }

    return ligands;
  };

  // Process the tree, from the bottom up
  for(auto it = byDepth.rbegin(); it != byDepth.rend(); ++it) {
    const auto& currentEdges = *it;

    for(const auto& edge: currentEdges) {
      auto sourceIndex = boost::source(edge, _tree);
      auto targetIndex = boost::target(edge, _tree);

      if(_doubleBondEdges.count(edge) > 0) {
        // Instantiate an EZStereocenter here!

        auto sourceIndicesToRank = _auxiliaryAdjacentsToRank(
          sourceIndex,
          {targetIndex}
        );

        auto targetIndicesToRank = _auxiliaryAdjacentsToRank(
          targetIndex,
          {sourceIndex}
        );

        // Ensure both index sets have either 1 or 2 elements
        if(
          !sourceIndicesToRank.empty() && sourceIndicesToRank.size() <= 2
          && !targetIndicesToRank.empty() && targetIndicesToRank.size() <= 2
        ) {
          auto makeEZRanking = [&](
            const TreeVertexIndex sourceIndex,
            const std::set<TreeVertexIndex>& indicesToRank
          ) -> RankingInformation {
            RankingInformation ranking;

            ranking.sortedSubstituents = _mapToAtomIndices(
              _auxiliaryApplySequenceRules(
                sourceIndex,
                indicesToRank
              )
            );

            ranking.ligands = auxiliaryLigands(ranking.sortedSubstituents);
            ranking.ligandsRanking = RankingInformation::rankLigands(
              ranking.ligands,
              ranking.sortedSubstituents
            );

            /* NOTE: There is no need to collect linking information since we
             * are in an acyclic digraph, and any cycles in the original molecule
             * are no longer present.
             */

            return ranking;
          };

          RankingInformation sourceRanking = makeEZRanking(
            sourceIndex,
            sourceIndicesToRank
          );

          RankingInformation targetRanking = makeEZRanking(
            targetIndex,
            targetIndicesToRank
          );

          const AtomIndexType molSourceIndex = _tree[sourceIndex].molIndex;
          const AtomIndexType molTargetIndex = _tree[targetIndex].molIndex;

          auto newStereocenter = Stereocenters::EZStereocenter {
            molSourceIndex,
            sourceRanking,
            molTargetIndex,
            targetRanking
          };

          // Try to assign the new stereocenter
          if(newStereocenter.numStereopermutations() > 1) {
            if(positionsOption) {
              // Fit from positions
              newStereocenter.fit(
                positionsOption.value()
              );
            } else {
              /* Need to get chiral information from the molecule (if present).
               * Have to be careful, any stereocenters on the same atom may have
               * different symmetries and different ranking for the same
               * substituents
               */
              if(_stereocentersRef.involving(molSourceIndex)) {
                const auto& stereocenterPtr = _stereocentersRef.at(molSourceIndex);

                if(
                  stereocenterPtr->type() == Stereocenters::Type::EZStereocenter
                  && stereocenterPtr->involvedAtoms() == std::vector<AtomIndexType> {
                    std::min(molSourceIndex, molTargetIndex),
                    std::max(molSourceIndex, molTargetIndex)
                  } && stereocenterPtr->numStereopermutations() == 2
                ) {
                  // Take the assignment from that stereocenter
                  newStereocenter.assign(
                    stereocenterPtr->assigned()
                  );
                }
              }
            }

            _tree[edge].stereocenterOption = newStereocenter;

            // Mark that we instantiated something
            foundEZStereocenters = true;
          }

          if /* C++17 constexpr */ (buildTypeIsDebug) {
            if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
              _writeGraphvizFiles({
                _adaptedMolGraphviz,
                dumpGraphviz("Sequence rule 3 preparation", {rootIndex})
              });
            }
          }
        }
      }

      if(
        !_tree[edge].stereocenterOption // No EZStereocenter on this edge
        && !_tree[targetIndex].stereocenterOption // No CNStereocenter
        && _nonDuplicateDegree(targetIndex) >= 3 // Min. degree for chirality
      ) {
        const AtomIndexType molSourceIndex = _tree[targetIndex].molIndex;

        RankingInformation centerRanking;

        centerRanking.sortedSubstituents = _mapToAtomIndices(
          _auxiliaryApplySequenceRules(
            targetIndex,
            _auxiliaryAdjacentsToRank(targetIndex, {})
          )
        );

        centerRanking.ligands = auxiliaryLigands(centerRanking.sortedSubstituents);
        centerRanking.ligandsRanking = RankingInformation::rankLigands(
          centerRanking.ligands,
          centerRanking.sortedSubstituents
        );

        // Again, no links since we're in an acyclic graph now
        Symmetry::Name localSymmetry;
        if(
          _stereocentersRef.involving(molSourceIndex)
          && _stereocentersRef.at(molSourceIndex)->type() == Stereocenters::Type::CNStereocenter
        ) {
          localSymmetry = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
            _stereocentersRef.at(molSourceIndex)
          ) -> getSymmetry();
        } else {
          localSymmetry = LocalGeometry::determineLocalGeometry(
            _graphRef,
            molSourceIndex,
            centerRanking
          );
        }
        const unsigned nHydrogens = _adjacentTerminalHydrogens(targetIndex);

        /* In case only one assignment is possible in the set of symmetries of
         * the same size, there is no reason to fit the stereocenter at all
         */
        if(
          temple::any_of(
            Symmetry::allNames,
            [&localSymmetry, nHydrogens](const Symmetry::Name name) -> bool {
              return (
                Symmetry::size(name) == Symmetry::size(localSymmetry)
                && Symmetry::hasMultipleUnlinkedAssignments(name, nHydrogens)
              );
            }
          )
        ) {
          // Instantiate a CNStereocenter here!

          auto newStereocenter = Stereocenters::CNStereocenter {
            _graphRef,
            localSymmetry,
            molSourceIndex,
            centerRanking
          };

          if(newStereocenter.numStereopermutations() > 1) {
            if(positionsOption) {
              pickyFit(
                newStereocenter,
                _graphRef,
                positionsOption.value(),
                localSymmetry
              );
            } else { // Try to get an assignment from the molecule
              if(_stereocentersRef.involving(molSourceIndex)) {
                const auto& stereocenterPtr = _stereocentersRef.at(molSourceIndex);

                /* TODO
                 * consider what to do in cases in which molecule number of
                 * assignments is fewer finding the equivalent case by means of
                 * rotations should be possible
                 *
                 * widening cases are not assignable
                 */
                if(
                  stereocenterPtr->type() == Stereocenters::Type::CNStereocenter
                  && stereocenterPtr->involvedAtoms() == std::vector<AtomIndexType> {
                    molSourceIndex
                  } && stereocenterPtr->numStereopermutations() == newStereocenter.numStereopermutations()
                ) {
                  newStereocenter.assign(
                    stereocenterPtr->assigned()
                  );
                }
              }
            }

            if(
              !disregardStereocenter(
                newStereocenter,
                graph::elementType(molSourceIndex, _graphRef),
                _cyclesRef,
                Options::temperatureRegime
              )
            ) {
              _tree[targetIndex].stereocenterOption = std::move(newStereocenter);

              // Mark that we instantiated something
              foundCNStereocenters = true;
            }
          }

          if /*C++17 constexpr */ (buildTypeIsDebug) {
            if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
              _writeGraphvizFiles({
                _adaptedMolGraphviz,
                dumpGraphviz("Sequence rule 3 preparation", {rootIndex})
              });
            }
          }
        }
      }
    }
  }

  // Was any EZStereocenter instantiated? If not, we can skip rule 3.
  if(foundEZStereocenters) {

    // Apply sequence rule 3
    _runBFS<
      3, // Sequence rule 3
      true, // BFS downwards only
      true, // Insert edges
      false, // Insert vertices
      TreeEdgeIndex, // Multiset value type
      SequenceRuleThreeEdgeComparator // Multiset comparator type
    >(
      rootIndex, // Source index is root
      _branchOrderingHelper
    );

    if /* C++17 constexpr */(buildTypeIsDebug) {
      Log::log(Log::Particulars::RankingTreeDebugInfo)
        << "Sets post sequence rule 3: {"
        << temple::condenseIterable(
          temple::map(
            _branchOrderingHelper.getSets(),
            [](const auto& indexSet) -> std::string {
              return "{"s + temple::condenseIterable(indexSet) + "}"s;
            }
          )
        ) << "}\n";
    }

    // Was sequence rule 3 enough?
    if(_branchOrderingHelper.isTotallyOrdered()) {
      return;
    }
  }

  // In case neither EZ or CNStereocenters were found, we can skip 4-5
  if(!foundEZStereocenters && !foundCNStereocenters) {
    return;
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
   *     Can probably use _auxiliaryApplySequenceRules recursively at the first
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

    std::map<
      TreeVertexIndex,
      std::multiset<VariantType, SequenceRuleFourVariantComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::vector<TreeVertexIndex>
    > seeds;

    undecidedSets = _branchOrderingHelper.getUndecidedSets();

    // Initialize the BFS state
    for(const auto& undecidedSet : undecidedSets) {
      for(const auto& undecidedBranch : undecidedSet) {
        comparisonSets.emplace(
          undecidedBranch,
          *this
        );

        comparisonSets.at(undecidedBranch).emplace(undecidedBranch);
        comparisonSets.at(undecidedBranch).emplace(
          boost::edge(rootIndex, undecidedBranch, _tree).first
        );

        seeds[undecidedBranch] = {undecidedBranch};
      }
    }

    /* First comparison (undecided sets are up-to-date from previous sequence
     * rule)
     */
    _compareBFSSets(comparisonSets, undecidedSets, _branchOrderingHelper);

    // Recalculate undecided sets
    undecidedSets = _branchOrderingHelper.getUndecidedSets();

    if /* C++17 constexpr */ (buildTypeIsDebug) {
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptedMolGraphviz,
          dumpGraphviz("Sequence rule 4A", {rootIndex}, _collectSeeds(seeds, undecidedSets)),
          _makeBFSStateGraph("4A", rootIndex, comparisonSets, undecidedSets),
          _branchOrderingHelper.dumpGraphviz()
        });
      }
    }

    while(!undecidedSets.empty() && _relevantSeeds(seeds, undecidedSets)) {
      // Perform a full BFS Step on all undecided set seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::vector<TreeVertexIndex> newSeeds;

          for(const auto& seed: branchSeeds) {
            for( // Out-edges
              auto outIterPair = boost::out_edges(seed, _tree);
              outIterPair.first != outIterPair.second;
              ++outIterPair.first
            ) {
              const auto& outEdge = *outIterPair.first;

              auto edgeTarget = boost::target(outEdge, _tree);

              if(!_tree[edgeTarget].isDuplicate) {
                comparisonSets.at(undecidedTreeIndex).emplace(edgeTarget);
                comparisonSets.at(undecidedTreeIndex).emplace(outEdge);

                newSeeds.push_back(edgeTarget);
              }
            }
          }

          // Overwrite the seeds
          branchSeeds = std::move(newSeeds);
        }
      }

      // Make comparisons in all undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, _branchOrderingHelper);

      // Recalculate the undecided sets
      undecidedSets = _branchOrderingHelper.getUndecidedSets();

      if /* C++17 constexpr */ (buildTypeIsDebug) {
        if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
          _writeGraphvizFiles({
            _adaptedMolGraphviz,
            dumpGraphviz("Sequence rule 4A", {rootIndex}, _collectSeeds(seeds, undecidedSets)),
            _makeBFSStateGraph("4A", rootIndex, comparisonSets, undecidedSets),
            _branchOrderingHelper.dumpGraphviz()
          });
        }
      }
    }

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "Sets post sequence rule 4A: {"
      << temple::condenseIterable(
        temple::map(
          _branchOrderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
    }

    // Is Sequence Rule 4, part A enough?
    if(_branchOrderingHelper.isTotallyOrdered()) {
      return;
    }

    /* Second part: Choosing representational stereodescriptors for all
     * branches and pairing in sequence of priority for branches.
     */
    std::map<
      TreeVertexIndex,
      std::set<VariantType>
    > stereocenterMap;

    VariantHasInstantiatedStereocenter isInstantiatedChecker {*this};

    // Copy all variants from part A that are actually instantiated
    for(const auto& undecidedSet : _branchOrderingHelper.getUndecidedSets()) {
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

      // If there are no stereocenters to rank, bounce
      if(stereocenterMap.at(branchIndex).empty()) {
        representativeStereodescriptors[branchIndex] = {};
        continue;
      }

      // Compare variants in a branch based on mixed depth
      temple::forAllPairs(
        variantSet,
        [&](const auto& a, const auto& b) {
          auto aDepth = boost::apply_visitor(depthFetcher, a);
          auto bDepth = boost::apply_visitor(depthFetcher, b);

          if(aDepth < bDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
          } else if(bDepth < aDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
          } else {
            // Try to resolve via ranking downwards from the junction
            auto junctionInfo = _junction(
              boost::apply_visitor(sourceNodeFetcher, a),
              boost::apply_visitor(sourceNodeFetcher, b)
            );

            auto aJunctionChild = junctionInfo.firstPath.back();
            auto bJunctionChild = junctionInfo.secondPath.back();

            /* Do not use _auxiliaryApplySequenceRules for root-level ranking,
             * it should only establish differences within branches
             */
            if(junctionInfo.junction != rootIndex) {
              auto relativeRank = _auxiliaryApplySequenceRules(
                junctionInfo.junction,
                {aJunctionChild, bJunctionChild},
                junctionInfo.firstPath.size() - 1
              );

              /* relativeRank can only have sizes 1 or 2, where size 1 means
               * that no difference was found
               */
              if(relativeRank.size() == 2) {
                if(relativeRank.front().front() == aJunctionChild) {
                  relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
                } else {
                  relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
                }
              }
            }
          }
        }
      );

      if /* C++17 constexpr */ (buildTypeIsDebug) {
        if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
          _writeGraphvizFiles({
            _adaptedMolGraphviz,
            dumpGraphviz("Sequence rule 4B prep", {rootIndex}),
            relativeOrders.at(branchIndex).dumpGraphviz()
          });
        }
      }

      // Pick the representative stereodescriptor for each branch!
      auto stereocenterSets = relativeOrders.at(branchIndex).getSets();
      if(stereocenterSets.back().size() == 1) {
        // 1 - The solitary highest ranked stereogenic group
        representativeStereodescriptors[branchIndex] = {
          stereocenterSets.back().front()
        };
      } else {
        // 2 - The stereodescriptor(s) occurring more often than all others
        auto groupedByStringRep = temple::groupByMapping(
          stereocenterSets.back(),
          [&](const auto& variantType) -> std::string {
            return boost::apply_visitor(stringRepFetcher, variantType);
          }
        );

        auto maxSize = temple::accumulate(
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
              *stringGroup.begin()
            );
          }
        }
      }
    }

    if /* C++17 constexpr */ (buildTypeIsDebug) {
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        std::set<TreeVertexIndex> representativeVertices;
        std::set<TreeEdgeIndex> representativeEdges;

        for(const auto& iterPair : representativeStereodescriptors) {
          const auto& representativeStereodescriptorsSet = iterPair.second;
          for(const auto& stereodescriptorVariant : representativeStereodescriptorsSet) {
            if(stereodescriptorVariant.which() == 0) {
              representativeVertices.insert(
                boost::get<TreeVertexIndex>(stereodescriptorVariant)
              );
            } else {
              representativeEdges.insert(
                boost::get<TreeEdgeIndex>(stereodescriptorVariant)
              );
            }
          }
        }

        _writeGraphvizFiles({
          _adaptedMolGraphviz,
          dumpGraphviz("Sequence rule 4B prep", {rootIndex}, representativeVertices, representativeEdges)
        });
      }
    }

    auto undecidedBranchSets = _branchOrderingHelper.getUndecidedSets();

    VariantLikePair variantLikeComparator {*this};

    for(const auto& undecidedSet : undecidedBranchSets) {
      temple::forAllPairs(
        undecidedSet,
        [&](const auto& branchA, const auto& branchB) {
          // Do nothing if neither have representative stereodescriptors
          if(
            representativeStereodescriptors.at(branchA).size() == 0
            && representativeStereodescriptors.at(branchB).size() == 0
          ) {
            return;
          }

          // Precedence via amount of representative stereodescriptors
          if(
            representativeStereodescriptors.at(branchA).size()
            < representativeStereodescriptors.at(branchB).size()
          ) {
            _branchOrderingHelper.addLessThanRelationship(branchA, branchB);
          } else if(
            representativeStereodescriptors.at(branchB).size()
            < representativeStereodescriptors.at(branchA).size()
          ) {
            _branchOrderingHelper.addLessThanRelationship(branchB, branchA);
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
             * This works because assignments are canonical, i.e. a
             * correspondence with R/S makes sense.
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
              temple::forAllPairs(
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
              temple::forAllPairs(
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
                _branchOrderingHelper.addLessThanRelationship(branchB, branchA);
                break;
              }

              if(BBranchLikePairs < ABranchLikePairs) {
                _branchOrderingHelper.addLessThanRelationship(branchA, branchB);
                break;
              }

              ++branchAStereocenterGroupIter;
              ++branchBStereocenterGroupIter;
            }

            if /* C++17 constexpr */ (buildTypeIsDebug) {
              if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
                _writeGraphvizFiles({
                  _adaptedMolGraphviz,
                  dumpGraphviz("4B"s),
                  _make4BGraph(
                    rootIndex,
                    representativeStereodescriptors,
                    branchA,
                    branchB,
                    branchAOrders,
                    branchBOrders,
                    branchAStereocenterGroupIter,
                    branchBStereocenterGroupIter
                  ),
                  _branchOrderingHelper.dumpGraphviz()
                });
              }
            }
          }
        }
      );
    }

    if /* C++17 constexpr */ (buildTypeIsDebug) {
      Log::log(Log::Particulars::RankingTreeDebugInfo)
        << "Sets post sequence rule 4B: {"
        << temple::condenseIterable(
          temple::map(
            _branchOrderingHelper.getSets(),
            [](const auto& indexSet) -> std::string {
              return "{"s + temple::condenseIterable(indexSet) + "}"s;
            }
          )
        ) << "}\n";
    }

    // Is Sequence Rule 4, part B enough?
    if(_branchOrderingHelper.isTotallyOrdered()) {
      return;
    }

    /* 4C Pseudo-asymmetries
     * This is skipped, since rule 5 covers this without issue
     */

  } // End sequence rule 4 local scope

  // Was sequence rule 4 enough?
  if(_branchOrderingHelper.isTotallyOrdered()) {
    return;
  }

  /* Sequence rule 5:
   * - Atom or group with {R, M, Z} precedes {S, P, E}
   */
  _runBFS<
    5, // Sequence rule 5
    true, // BFS downwards only
    true, // Insert edges
    true, // Insert vertices
    VariantType, // MultisetValueType
    SequenceRuleFiveVariantComparator // Multiset comparator type
  >(
    rootIndex, // Source index is root
    _branchOrderingHelper
  );

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "Sets post sequence rule 5: {"
      << temple::condenseIterable(
        temple::map(
          _branchOrderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
  }

  // Exhausted sequence rules, anything undecided is now equal
}

std::vector<
  std::vector<RankingTree::TreeVertexIndex>
> RankingTree::_auxiliaryApplySequenceRules(
  const RankingTree::TreeVertexIndex& sourceIndex,
  const std::set<RankingTree::TreeVertexIndex>& adjacentsToRank,
  const boost::optional<unsigned>& depthLimitOptional
) const {
  /* Sequence rule 1 */
  OrderDiscoveryHelper<TreeVertexIndex> orderingHelper {adjacentsToRank};

  // Return immediately if depthLimit is zero
  if(depthLimitOptional && depthLimitOptional.value() == 0u) {
    return orderingHelper.getSets();
  }

#ifdef RANKING_TREE_OPTIMIZATION_REUSE_AUXILIARY_RESULTS
  orderingHelper.addRelationshipsFromOther(_allOrdering);
  // Is the information from _allOrdering sufficient?
  if(orderingHelper.isTotallyOrdered()) {
    return orderingHelper.getSets();
  }
#endif

#ifndef NEDEBUG
  Log::log(Log::Particulars::RankingTreeDebugInfo)
    << "  Auxiliary ranking substituents of tree index "
    << sourceIndex
    <<  ": "
    << temple::condenseIterable(
      temple::map(
        orderingHelper.getSets(),
        [](const auto& indexSet) -> std::string {
          return "{"s + temple::condenseIterable(indexSet) + "}"s;
        }
      )
    ) << "\n";
#endif


  auto undecidedSets = orderingHelper.getUndecidedSets();

  _runBFS<
    1, // Sequence rule 1
    false, // BFS downwards only
    false, // Insert edges
    true, // Insert vertices
    TreeVertexIndex, // Multiset value type
    SequenceRuleOneVertexComparator // Multiset comparator type
  >(
    sourceIndex,
    orderingHelper,
    depthLimitOptional
  );

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "  Sets post sequence rule 1: {"
      << temple::condenseIterable(
        temple::map(
          orderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
  }

  // Is Sequence Rule 1 enough?
  if(orderingHelper.isTotallyOrdered()) {
#ifdef RANKING_TREE_OPTIMIZATION_REUSE_AUXILIARY_RESULTS
    if(!depthLimitOptional) {
      _allOrdering.addAllFromOther(orderingHelper);
    }
#endif
    // No conversion of indices in _auxiliaryApplySequenceRules()!
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
  _runBFS<
    3, // Sequence rule 3
    false, // BFS downwards only
    true, // Insert edges
    false, // Insert vertices
    TreeEdgeIndex, // Multiset value type
    SequenceRuleThreeEdgeComparator // Multiset comparator type
  >(
    sourceIndex,
    orderingHelper,
    depthLimitOptional
  );

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "  Sets post sequence rule 3: {"
      << temple::condenseIterable(
        temple::map(
          orderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
  }

  // Is sequence rule 3 enough?
  if(orderingHelper.isTotallyOrdered()) {
#ifdef RANKING_TREE_OPTIMIZATION_REUSE_AUXILIARY_RESULTS
    if(!depthLimitOptional) {
      _allOrdering.addAllFromOther(orderingHelper);
    }
#endif
    // No conversion of indices in _auxiliaryApplySequenceRules()!
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
   *     Can probably use _auxiliaryApplySequenceRules recursively at the first
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

    std::set<TreeVertexIndex> visitedVertices {sourceIndex};

    std::map<
      TreeVertexIndex,
      std::multiset<VariantType, SequenceRuleFourVariantComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::vector<TreeVertexIndex>
    > seeds;

    undecidedSets = orderingHelper.getUndecidedSets();

    // Initialize the BFS state
    for(const auto& undecidedSet : undecidedSets) {
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

    _compareBFSSets(comparisonSets, undecidedSets, orderingHelper);

    undecidedSets = orderingHelper.getUndecidedSets();

    unsigned depth = 1;

    if /* C++17 constexpr */ (buildTypeIsDebug) {
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        _writeGraphvizFiles({
          _adaptedMolGraphviz,
          dumpGraphviz("aux Sequence rule 4A", {sourceIndex}, visitedVertices),
          _makeBFSStateGraph("aux 4A", sourceIndex, comparisonSets, undecidedSets),
          orderingHelper.dumpGraphviz()
        });
      }
    }

    while(
      !undecidedSets.empty()
      && _relevantSeeds(seeds, undecidedSets)
      && depth < depthLimitOptional.value_or(std::numeric_limits<unsigned>::max())
    ) {
      // Perform a full BFS Step on all undecided set seeds
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedTreeIndex : undecidedSet) {
          auto& branchSeeds = seeds[undecidedTreeIndex];

          std::vector<TreeVertexIndex> newSeeds;

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

                newSeeds.push_back(edgeSource);
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

                newSeeds.push_back(edgeTarget);
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

      // Increment depth
      ++depth;

      if /* C++17 constexpr */ (buildTypeIsDebug) {
        if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
          _writeGraphvizFiles({
            _adaptedMolGraphviz,
            dumpGraphviz("aux Sequence rule 4A", {sourceIndex}, visitedVertices),
            _makeBFSStateGraph("aux 4A", sourceIndex, comparisonSets, undecidedSets),
            orderingHelper.dumpGraphviz()
          });
        }
      }
    }

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "  Sets post sequence rule 4A: {"
      << temple::condenseIterable(
        temple::map(
          orderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
  }

    // Is Sequence Rule 4, part A enough?
    if(orderingHelper.isTotallyOrdered()) {
#ifdef RANKING_TREE_OPTIMIZATION_REUSE_AUXILIARY_RESULTS
      if(!depthLimitOptional) {
        _allOrdering.addAllFromOther(orderingHelper);
      }
#endif
      // No conversion of indices in _auxiliaryApplySequenceRules()!
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

      // If there are no stereocenters to rank, bounce
      if(stereocenterMap.at(branchIndex).empty()) {
        representativeStereodescriptors[branchIndex] = {};
        continue;
      }

      // Compare based on depth
      temple::forAllPairs(
        variantSet,
        [&](const auto& a, const auto& b) {
          auto aDepth = boost::apply_visitor(depthFetcher, a);
          auto bDepth = boost::apply_visitor(depthFetcher, b);

          if(aDepth < bDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
          } else if(bDepth < aDepth) {
            relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
          } else {
            auto junctionInfo = _junction(
              boost::apply_visitor(sourceNodeFetcher, a),
              boost::apply_visitor(sourceNodeFetcher, b)
            );

            auto aJunctionChild = junctionInfo.firstPath.back();
            auto bJunctionChild = junctionInfo.secondPath.back();

            /* Do not use _auxiliaryApplySequenceRules for root-level ranking,
             * it should only establish differences within branches here
             */
            if(junctionInfo.junction != rootIndex) {
              auto relativeRank = _auxiliaryApplySequenceRules(
                junctionInfo.junction,
                {aJunctionChild, bJunctionChild},
                junctionInfo.firstPath.size() - 1
              );

              /* relativeRank can only have sizes 1 or 2, 1 meaning that no
               * difference was found
               */
              if(relativeRank.size() == 2) {
                if(relativeRank.front().front() == aJunctionChild) {
                  relativeOrders.at(branchIndex).addLessThanRelationship(b, a);
                } else {
                  relativeOrders.at(branchIndex).addLessThanRelationship(a, b);
                }
              }
            }
          }
        }
      );

      if /* C++17 constexpr */ (buildTypeIsDebug) {
        if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
          _writeGraphvizFiles({
            _adaptedMolGraphviz,
            dumpGraphviz("Aux 4B prep", {rootIndex}),
            relativeOrders.at(branchIndex).dumpGraphviz()
          });
        }
      }

      // Pick the representative stereodescriptor for each branch!
      auto stereocenterSets = relativeOrders.at(branchIndex).getSets();

      if(stereocenterSets.back().size() == 1) {
        // 1 - The solitary highest ranked stereogenic group
        representativeStereodescriptors[branchIndex] = {
          stereocenterSets.back().front()
        };
      } else {
        // 2 - The stereodescriptor occurring more often than all others
        auto groupedByStringRep = temple::groupByMapping(
          stereocenterSets.back(),
          [&](const auto& variantType) -> std::string {
            return boost::apply_visitor(stringRepFetcher, variantType);
          }
        );

        auto maxSize = temple::accumulate(
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
              *stringGroup.begin()
            );
          }
        }
      }
    }

    auto undecidedBranchSets = orderingHelper.getUndecidedSets();

    VariantLikePair variantLikeComparator {*this};

    for(const auto& undecidedSet : undecidedBranchSets) {
      temple::forAllPairs(
        undecidedSet,
        [&](const auto& branchA, const auto& branchB) {
          // Do nothing if neither have representative stereodescriptors
          if(
            representativeStereodescriptors.at(branchA).size() == 0
            && representativeStereodescriptors.at(branchB).size() == 0
          ) {
            return;
          }

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
             * This works because assignments are canonical, i.e. a
             * correspondence with R/S makes sense.
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
              temple::forAllPairs(
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
              temple::forAllPairs(
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
                orderingHelper.addLessThanRelationship(branchB, branchA);
                break;
              }

              if(BBranchLikePairs < ABranchLikePairs) {
                orderingHelper.addLessThanRelationship(branchA, branchB);
                break;
              }

              ++branchAStereocenterGroupIter;
              ++branchBStereocenterGroupIter;
            }

            if /* C++17 constexpr */ (buildTypeIsDebug) {
              if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
                _writeGraphvizFiles({
                  _adaptedMolGraphviz,
                  dumpGraphviz("aux 4B"s),
                  _make4BGraph(
                    sourceIndex,
                    representativeStereodescriptors,
                    branchA,
                    branchB,
                    branchAOrders,
                    branchBOrders,
                    branchAStereocenterGroupIter,
                    branchBStereocenterGroupIter
                  ),
                  orderingHelper.dumpGraphviz()
                });
              }
            }
          }
        }
      );
    }

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "  Sets post sequence rule 4B: {"
      << temple::condenseIterable(
        temple::map(
          orderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
  }

    // Is Sequence Rule 4, part B enough?
    if(orderingHelper.isTotallyOrdered()) {
#ifdef RANKING_TREE_OPTIMIZATION_REUSE_AUXILIARY_RESULTS
      if(!depthLimitOptional) {
        _allOrdering.addAllFromOther(orderingHelper);
      }
#endif
      // No conversion of indices in _auxiliaryApplySequenceRules()!
      return orderingHelper.getSets();
    }

    /* 4C Pseudo-asymmetries
     * This is skipped, since rule 5 covers this without issue
     */

  } // End sequence rule 4 local scope

  /* Sequence rule 5:
   * - Atom or group with {R, M, Z} precedes {S, P, E}
   */
  _runBFS<
    5, // Sequence rule 5
    false, // BFS downwards only
    true, // Insert edges
    true, // Insert vertices
    VariantType, // MultisetValueType
    SequenceRuleFiveVariantComparator // Multiset comparator type
  >(
    sourceIndex,
    orderingHelper,
    depthLimitOptional
  );

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "  Sets post sequence rule 5: {"
      << temple::condenseIterable(
        temple::map(
          orderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
  }

#ifdef RANKING_TREE_OPTIMIZATION_REUSE_AUXILIARY_RESULTS
  if(!depthLimitOptional) {
    _allOrdering.addAllFromOther(orderingHelper);
  }
#endif
  // Exhausted sequence rules, return the sets
  return orderingHelper.getSets();
}

std::vector<RankingTree::TreeVertexIndex> RankingTree::_expand(
  const RankingTree::TreeVertexIndex& index,
  const std::set<AtomIndexType>& molIndicesInBranch
) {
  std::vector<TreeVertexIndex> newIndices;

  std::set<AtomIndexType> treeOutAdjacencies;
  TreeGraphType::out_edge_iterator iter, end;
  std::tie(iter, end) = boost::out_edges(index, _tree);
  while(iter != end) {
    auto targetVertex = boost::target(*iter, _tree);

    treeOutAdjacencies.insert(
      _tree[targetVertex].molIndex
    );

    newIndices.push_back(targetVertex);

    ++iter;
  }

  for(
    const auto& molAdjacentIndex :
    graph::adjacents(_tree[index].molIndex, _graphRef)
  ) {
    if(treeOutAdjacencies.count(molAdjacentIndex) != 0) {
      continue;
    }

    if(molIndicesInBranch.count(molAdjacentIndex) == 0) {
      auto newIndex = boost::add_vertex(_tree);
      _tree[newIndex].molIndex = molAdjacentIndex;
      _tree[newIndex].isDuplicate = false;

      newIndices.push_back(newIndex);

      boost::add_edge(index, newIndex, _tree);

      // Need to add duplicates!
      auto duplicateVertices = _addBondOrderDuplicates(index, newIndex);
      std::copy(
        duplicateVertices.begin(),
        duplicateVertices.end(),
        std::back_inserter(newIndices)
      );

      // Copy the indices and insert the newest addition
      auto molIndicesInBranchCopy = molIndicesInBranch;
      molIndicesInBranchCopy.insert(molAdjacentIndex);
    } else if(molAdjacentIndex != _tree[_parent(index)].molIndex)  {
      auto newIndex = boost::add_vertex(_tree);
      _tree[newIndex].molIndex = molAdjacentIndex;
      _tree[newIndex].isDuplicate = true;

      newIndices.push_back(newIndex);

      boost::add_edge(index, newIndex, _tree);
    }
  }

  // Does not contain any new duplicate indices
  return newIndices;
}

//! Specialization of toString for tree vertices
template<>
std::string RankingTree::toString(const TreeVertexIndex& vertexIndex) const {
  return std::to_string(vertexIndex);
}

//! Specialization of toString for tree edges
template<>
std::string RankingTree::toString(const TreeEdgeIndex& edgeIndex) const {
  return (
    std::to_string(
      boost::source(edgeIndex, _tree)
    ) + "→"s + std::to_string(
      boost::target(edgeIndex, _tree)
    )
  );
}

//! Specialization of toString for variant types containing either a vertex or edge
template<>
std::string RankingTree::toString(const VariantType& variant) const {
  if(variant.which() == 0) {
    return toString<TreeVertexIndex>(boost::get<TreeVertexIndex>(variant));
  }

  return toString<TreeEdgeIndex>(boost::get<TreeEdgeIndex>(variant));
}

//! Returns the parent of a node. Fails if called on the root!
RankingTree::TreeVertexIndex RankingTree::_parent(const RankingTree::TreeVertexIndex& index) const {
  assert(index != rootIndex);

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

unsigned RankingTree::_adjacentTerminalHydrogens(const TreeVertexIndex& index) const {
  unsigned count = 0;

  auto outIterPair = boost::out_edges(index, _tree);
  while(outIterPair.first != outIterPair.second) {
    auto edgeTarget = boost::target(*outIterPair.first, _tree);

    if(
      graph::elementType(_tree[edgeTarget].molIndex, _graphRef) == Delib::ElementType::H
      && boost::out_degree(edgeTarget, _tree) == 0
    ) {
      ++count;
    }

    ++outIterPair.first;
  }

  return count;
}

bool RankingTree::_isBondSplitDuplicateVertex(const TreeVertexIndex& index) const {
  // If the vertex isn't duplicate, it can't be a closure
  if(!_tree[index].isDuplicate) {
    return false;
  }

  /* Duplicate atoms have two sources - they can be either from multiple bond
   * breakage or from cycle closure. If we can find the typical pattern from a
   * multiple bond breakage, then it can't be a cycle closure
   */
  auto parentIndex = _parent(index);

  auto inEdgesIterPair = boost::in_edges(parentIndex, _tree);
  while(inEdgesIterPair.first != inEdgesIterPair.second) {
    auto adjacentIndex = boost::source(*inEdgesIterPair.first, _tree);
    if(
      _tree[adjacentIndex].molIndex == _tree[index].molIndex
      && !_tree[adjacentIndex].isDuplicate
    ) {
      return true;
    }

    ++inEdgesIterPair.first;
  }

  auto outEdgesIterPair = boost::out_edges(parentIndex, _tree);
  while(outEdgesIterPair.first != outEdgesIterPair.second) {
    auto adjacentIndex = boost::target(*outEdgesIterPair.first, _tree);

    if(
      _tree[adjacentIndex].molIndex == _tree[index].molIndex
      && !_tree[adjacentIndex].isDuplicate
    ) {
      return true;
    }

    ++outEdgesIterPair.first;
  }

  // Remaining case is a cycle closure
  return false;
}

bool RankingTree::_isCycleClosureDuplicateVertex(const TreeVertexIndex& index) const {
  // If the vertex isn't duplicate, it can't be a closure
  if(!_tree[index].isDuplicate) {
    return false;
  }

  /* Duplicate atoms have two sources - they can be either from multiple bond
   * breakage or from cycle closure. If we can find the typical pattern from a
   * multiple bond breakage, then it isn't a cycle closure.
   */
  auto parentIndex = _parent(index);

  auto inEdgesIterPair = boost::in_edges(parentIndex, _tree);
  while(inEdgesIterPair.first != inEdgesIterPair.second) {
    auto adjacentIndex = boost::source(*inEdgesIterPair.first, _tree);
    if(
      _tree[adjacentIndex].molIndex == _tree[index].molIndex
      && !_tree[adjacentIndex].isDuplicate
    ) {
      return false;
    }

    ++inEdgesIterPair.first;
  }

  auto outEdgesIterPair = boost::out_edges(parentIndex, _tree);
  while(outEdgesIterPair.first != outEdgesIterPair.second) {
    auto adjacentIndex = boost::target(*outEdgesIterPair.first, _tree);

    if(
      _tree[adjacentIndex].molIndex == _tree[index].molIndex
      && !_tree[adjacentIndex].isDuplicate
    ) {
      return false;
    }

    ++outEdgesIterPair.first;
  }

  /* When you have eliminated the impossible, whatever remains, however
   * improbable, must be the truth.
   */
  return true;
}

std::set<RankingTree::TreeVertexIndex> RankingTree::_auxiliaryAdjacentsToRank(
  const TreeVertexIndex& sourceIndex,
  const std::set<TreeVertexIndex>& excludeIndices
) const {
  return temple::moveIf(
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
         * to be kept, while multiple bond splits need to be removed.
         */
        return _isCycleClosureDuplicateVertex(nodeIndex);
      }

      // Keep the vertex in all other cases
      return true;
    }
  );
}

unsigned RankingTree::_nonDuplicateDegree(const RankingTree::TreeVertexIndex& index) const {
  auto adjacents = _adjacents(index);

  auto numDuplicate = temple::accumulate(
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
    std::vector<RankingTree::TreeVertexIndex>
  >& seeds,
  const std::vector<
    std::vector<RankingTree::TreeVertexIndex>
  >& undecidedSets
) {
  for(const auto& undecidedSet : undecidedSets) {
    for(const auto& undecidedTreeIndex : undecidedSet) {
      if(!seeds.at(undecidedTreeIndex).empty()) {
        return true;
      }
    }
  }

  return false;
}

std::vector<RankingTree::TreeVertexIndex> RankingTree::_addBondOrderDuplicates(
  const TreeVertexIndex& treeSource,
  const TreeVertexIndex& treeTarget
) {
  std::vector<TreeVertexIndex> newIndices;

  auto molGraphEdge = graph::edge(
    _tree[treeSource].molIndex,
    _tree[treeTarget].molIndex,
    _graphRef
  );

  /* In case the bond order is non-fractional (aromatic / eta)
   * and > 1, add duplicate atoms
   */

  auto bondType = graph::bondType(molGraphEdge, _graphRef);

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
  auto integralBondOrder = static_cast<unsigned>(bondOrder);

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

      newIndices.push_back(b);

      boost::add_edge(treeSource, b, _tree);
    }
  }

  return newIndices;
}

std::set<RankingTree::TreeVertexIndex> RankingTree::_treeIndicesInBranch(TreeVertexIndex index) const {
  std::set<AtomIndexType> indices {index};

  while(index != rootIndex) {
    // All nodes in the graph must have in_degree of 1
    assert(boost::in_degree(index, _tree) == 1);

    // Follow singular in_node to source, overwrite index
    auto iterPair = boost::in_edges(index, _tree);
    index = boost::source(*iterPair.first, _tree);

    indices.insert(index);
  }

  return indices;
}

std::set<AtomIndexType> RankingTree::_molIndicesInBranch(const TreeVertexIndex& index) const {
  return temple::map(
    _treeIndicesInBranch(index),
    [&](const auto& treeIndex) -> AtomIndexType {
      return _tree[treeIndex].molIndex;
    }
  );
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

  while(index != rootIndex) {
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

  while(index != rootIndex) {
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

//! A data class, stores both the junction vertex and paths from the source vertices.
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
  std::set<TreeVertexIndex> aBranchIndices = _treeIndicesInBranch(a);

  // By default, the junction is root
  data.junction = rootIndex;

  // In case b is included in the path, that is the junction
  if(aBranchIndices.count(b) > 0) {
    data.junction = b;
  } else {
    // Backtrack with b, checking at each iteration
    auto bCurrent = b;
    while(bCurrent != rootIndex) {
      bCurrent = _parent(bCurrent);

      if(aBranchIndices.count(bCurrent) > 0) {
        data.junction = bCurrent;
        break;
      }
    }
  }

  /* Reconstruct the paths to the junction */

  for(
    TreeVertexIndex aCurrent = a;
    aCurrent != data.junction;
    aCurrent = _parent(aCurrent)
  ) {
    data.firstPath.push_back(aCurrent);
  }

  for(
    TreeVertexIndex bCurrent = b;
    bCurrent != data.junction;
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

  while(treeIndex != rootIndex) {
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


std::set<RankingTree::TreeVertexIndex> RankingTree::_collectSeeds(
  const std::map<
    TreeVertexIndex,
    std::vector<TreeVertexIndex>
  >& seeds,
  const std::vector<
    std::vector<TreeVertexIndex>
  >& undecidedSets
) {
  std::set<TreeVertexIndex> visited;

  for(const auto& undecidedSet : undecidedSets) {
    for(const auto& undecidedBranch : undecidedSet) {
      for(const auto& seed : seeds.at(undecidedBranch)) {
        visited.insert(seed);
      }
    }
  }

  return visited;
}

RankingTree::RankingTree(
  const GraphType& graph,
  const Cycles& cycles,
  const StereocenterList& stereocenters,
  std::string molGraphviz,
  const AtomIndexType atomToRank,
  const std::set<AtomIndexType>& excludeIndices,
  const ExpansionOption& expansionMethod,
  const boost::optional<AngstromWrapper>& positionsOption
) : _graphRef(graph),
    _cyclesRef(cycles),
    _stereocentersRef(stereocenters),
    _adaptedMolGraphviz(_adaptMolGraph(std::move(molGraphviz)))
{
  // Set the root vertex
  auto rootIndex = boost::add_vertex(_tree);
  _tree[rootIndex].molIndex = atomToRank;
  _tree[rootIndex].isDuplicate = false;

  /* Add the direct descendants of the root atom, if they are not explicitly
   * excluded by parameters
   */
  std::set<TreeVertexIndex> branchIndices;
  for(const auto rootAdjacentIndex : graph::adjacents(atomToRank, _graphRef)) {
    if(excludeIndices.count(rootAdjacentIndex) == 0) {
      auto branchIndex = boost::add_vertex(_tree);
      _tree[branchIndex].molIndex = rootAdjacentIndex;
      _tree[branchIndex].isDuplicate = false;

      boost::add_edge(rootIndex, branchIndex, _tree);

      _addBondOrderDuplicates(rootIndex, branchIndex);

      branchIndices.insert(branchIndex);
    }
  }

  // Set the ordering helper's list of unordered values
  _branchOrderingHelper.setUnorderedValues(branchIndices);

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "Ranking substituents of atom index "
      << _tree[rootIndex].molIndex
      << ": "
      << temple::condenseIterable(
        temple::map(
          _branchOrderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "\n";
  }

  // The class should work regardless of which tree expansion method is chosen
  if(expansionMethod == ExpansionOption::Optimized) {
    /* Expand only those branches which are yet undecided, immediately compare
     * using sequence rule 1.
     */

    std::map<
      TreeVertexIndex,
      std::multiset<TreeVertexIndex, SequenceRuleOneVertexComparator>
    > comparisonSets;

    std::map<
      TreeVertexIndex,
      std::vector<TreeVertexIndex>
    > seeds;

    for(const auto& branchIndex : branchIndices) {
      comparisonSets.emplace(
        branchIndex,
        *this
      );

      comparisonSets.at(branchIndex).insert(branchIndex);

      seeds[branchIndex] = {branchIndex};
    }

    auto undecidedSets = _branchOrderingHelper.getUndecidedSets();

    // Make the first comparison
    _compareBFSSets(comparisonSets, undecidedSets, _branchOrderingHelper);

    // Update the undecided sets
    undecidedSets = _branchOrderingHelper.getUndecidedSets();

    if /* C++17 constexpr */ (buildTypeIsDebug) {
      // Write debug graph files if the corresponding log particular is set
      if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
        std::string header = "Sequence rule 1";

        _writeGraphvizFiles({
          _adaptedMolGraphviz,
          dumpGraphviz(
            header,
            {rootIndex},
            _collectSeeds(seeds, undecidedSets)
          ),
          _makeBFSStateGraph(
            "R1"s,
            rootIndex,
            comparisonSets,
            undecidedSets
          ),
          _branchOrderingHelper.dumpGraphviz()
        });
      }
    }

    // TODO try to avoid repeated computation with _molIndicesInBranch somehow
    // Main BFS loop
    while(!undecidedSets.empty() && _relevantSeeds(seeds, undecidedSets)) {

      // Perform a step
      for(const auto& undecidedSet : undecidedSets) {
        for(const auto& undecidedBranch : undecidedSet) {
          auto& branchSeeds = seeds.at(undecidedBranch);
          auto& branchComparisonSet = comparisonSets.at(undecidedBranch);

          branchComparisonSet.clear();

          std::vector<TreeVertexIndex> newSeeds;

          for(const auto& seedVertex : branchSeeds) {
            auto newVertices = _expand(
              seedVertex,
              _molIndicesInBranch(seedVertex)
            );

            // Add ALL new vertices to the comparison set
            std::copy(
              newVertices.begin(),
              newVertices.end(),
              std::inserter(branchComparisonSet, branchComparisonSet.begin())
            );

            // Add non-duplicate vertices to the new seeds
            std::copy_if(
              newVertices.begin(),
              newVertices.end(),
              std::back_inserter(newSeeds),
              [&](const auto& newVertex) -> bool {
                return !_tree[newVertex].isDuplicate;
              }
            );
          }

          branchSeeds = std::move(newSeeds);
        }
      }

      // Compare and update the undecided sets
      _compareBFSSets(comparisonSets, undecidedSets, _branchOrderingHelper);
      undecidedSets = _branchOrderingHelper.getUndecidedSets();

      if /* C++17 constexpr */ (buildTypeIsDebug) {
        // Write debug graph files if the corresponding log particular is set
        if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
          std::string header = "Sequence rule 1";

          _writeGraphvizFiles({
            _adaptedMolGraphviz,
            dumpGraphviz(
              header,
              {rootIndex},
              _collectSeeds(seeds, undecidedSets)
            ),
            _makeBFSStateGraph(
              "R1"s,
              rootIndex,
              comparisonSets,
              undecidedSets
            ),
            _branchOrderingHelper.dumpGraphviz()
          });
        }
      }
    }

  } else { // Full tree expansion requested

    std::vector<TreeVertexIndex> seeds;
    std::copy(
      branchIndices.begin(),
      branchIndices.end(),
      std::back_inserter(seeds)
    );

    do {
      std::vector<TreeVertexIndex> newSeeds;

      for(const auto& seedVertex : seeds) {
        auto newVertices = _expand(
          seedVertex,
          _molIndicesInBranch(seedVertex)
        );

        // Add non-duplicate vertices to the new seeds
        std::copy_if(
          newVertices.begin(),
          newVertices.end(),
          std::back_inserter(newSeeds),
          [&](const auto& newVertex) -> bool {
            return !_tree[newVertex].isDuplicate;
          }
        );
      }

      seeds = std::move(newSeeds);
    } while(!seeds.empty());

    /* Apply sequence rule 1 here, not in _applySequenceRules, since the
     * optimized version IS the application of sequence rule 1
     */
    /* Sequence rule 1
     * - Higher atomic number precedes lower atomic number
     * - A duplicate atom node whose corresponding non-duplicate atom node is
     *   root or is closer to the root preceds a duplicate node whose
     *   corresponding atom node is further from the root
     */
    _runBFS<
      1, // Sequence rule 1
      true, // BFS downwards only
      false, // Insert edges
      true, // Insert vertices
      TreeVertexIndex, // Multiset value type
      SequenceRuleOneVertexComparator // Multiset comparator type
    >(
      rootIndex, // Source index is root
      _branchOrderingHelper
    );
  }

  if /* C++17 constexpr */ (buildTypeIsDebug) {
    Log::log(Log::Particulars::RankingTreeDebugInfo)
      << "Sets post sequence rule 1: {"
      << temple::condenseIterable(
        temple::map(
          _branchOrderingHelper.getSets(),
          [](const auto& indexSet) -> std::string {
            return "{"s + temple::condenseIterable(indexSet) + "}"s;
          }
        )
      ) << "}\n";
  }

  // Was sequence rule 1 enough?
  if(_branchOrderingHelper.isTotallyOrdered()) {
    return;
  }

  // Perform ranking
  _applySequenceRules(positionsOption);

#ifdef RANKING_TREE_OPTIMIZATION_REUSE_AUXILIARY_RESULTS
  if /* C++17 constexpr */ (buildTypeIsDebug) {
    if(Log::particulars.count(Log::Particulars::RankingTreeDebugInfo) > 0) {
      _writeGraphvizFiles({
        _adaptMolGraph(_moleculeRef.dumpGraphviz()),
        dumpGraphviz(
          "Final"s,
          {rootIndex}
        ),
        _branchOrderingHelper.dumpGraphviz(),
        _allOrdering.dumpGraphviz()
      });
    }
  }
#endif
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

std::string RankingTree::_make4BGraph(
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
) const {
  VariantStereocenterStringRepresentation stringFetcher {*this};
  VariantLikePair likeComparator {*this};

  const std::string nl = "\n"s;
  const std::string tableBegin = R"(<<table border="0" cellspacing="0" cellpadding="4">)";
  const std::string tableEnd = R"(</table>>)";
  const std::string rowBegin = R"(<tr>)";
  const std::string rowEnd = R"(</tr>)";
  const std::string br = R"(<br />)";
  const std::string cellEnd = R"(</td>)";
  const std::string greenColor = "forestgreen";
  const std::string redColor = "orangered";
  const std::string branchColor = "gray60";

  auto cellBegin = [&](
    const unsigned& span = 1,
    const boost::optional<std::string>& colorOption = boost::none
  ) -> std::string {
    std::string html = R"(<td border="1")";

    if(span != 1) {
      html += R"( colspan=")" + std::to_string(span) + R"(")";
    }

    if(colorOption) {
      html += R"( bgcolor=")" + colorOption.value() + R"(")";
    }

    html += ">"s;

    return html;
  };

  std::string graphviz = (
    "digraph G {\n"s
    + R"(  graph [fontname="Arial"];)" + nl
    + R"(  node [fontname="Arial", shape="record"];)" + nl
    + R"(  edge [fontname="Arial"];)" + nl
  );

  // Root node
  graphviz += R"(  root [)";
  graphviz += R"(label="4B\n\n)" + std::to_string(sourceIndex) + R"(")";
  graphviz += R"(, shape="square")";
  graphviz += R"(];)" + nl;

  auto makeBranchNode = [&](
    const auto& nodeName,
    const auto& branchIndex
  ) {
    std::string nodeGraphviz = "  "s + nodeName + "["s;
    nodeGraphviz += R"(shape="none", label=)" + tableBegin + rowBegin;
    for(const auto& representativeVariant : representativeStereodescriptors.at(branchIndex)) {
      nodeGraphviz += cellBegin(
        representativeStereodescriptors.at(branchIndex).size(),
        branchColor
      );
      nodeGraphviz += toString(representativeVariant) + br;
      nodeGraphviz += boost::apply_visitor(stringFetcher, representativeVariant) + cellEnd;
    }
    nodeGraphviz += rowEnd + tableEnd + R"(];)" + nl;

    return nodeGraphviz;
  };

  graphviz += makeBranchNode("branchA", branchA);
  graphviz += makeBranchNode("branchB", branchB);

  graphviz += R"(  root -> branchA;)" + nl;
  graphviz += R"(  root -> branchB;)" + nl;

  // Branch nodes
  auto makeBranch = [&](
    auto lastNode,
    const auto& nodePrefix,
    const auto& branchIndex,
    const auto& branchOrders,
    const auto& branchIter
  ) {
    std::string branchGraphviz;
    unsigned i = 0;

    for(auto iter = branchOrders.rbegin(); iter != branchOrders.rend(); ++iter) {
      const auto& variantList = *iter;

      // Node
      branchGraphviz += "  "s + nodePrefix + std::to_string(i) + " [";
      branchGraphviz += R"(shape="none", label=)" + tableBegin;

      std::string stereocenterRow = rowBegin;
      std::string likeRow = rowBegin;

      for(const auto& stereocenterVariant : variantList) {
        stereocenterRow += cellBegin(
          representativeStereodescriptors.at(branchIndex).size()
        );
        stereocenterRow += toString(stereocenterVariant) + br;
        stereocenterRow += boost::apply_visitor(stringFetcher, stereocenterVariant);
        stereocenterRow += cellEnd;

        if(branchIter == iter) {
          for(const auto& representativeVariant : representativeStereodescriptors.at(branchIndex)) {
            if(boost::apply_visitor(likeComparator, stereocenterVariant, representativeVariant)) {
              likeRow += cellBegin(1, greenColor) + cellEnd;
            } else {
              likeRow += cellBegin(1, redColor) + cellEnd;
            }
          }
        }
      }

      branchGraphviz += stereocenterRow + rowEnd;

      if(branchIter == iter) {
        branchGraphviz += likeRow + rowEnd;
      }

      branchGraphviz += tableEnd + R"(];)";
      branchGraphviz += nl;

      // Edge
      branchGraphviz += "  "s + lastNode + " -> " + nodePrefix + std::to_string(i) + ";"s + nl;

      lastNode = nodePrefix + std::to_string(i);
      ++i;
    }

    return branchGraphviz;
  };

  graphviz += makeBranch("branchA"s, "a"s, branchA, branchAOrders, branchAIter);
  graphviz += makeBranch("branchB"s, "b"s, branchB, branchBOrders, branchBIter);

  graphviz += "}";

  return graphviz;
}

std::vector<
  std::vector<AtomIndexType>
> RankingTree::_mapToAtomIndices(
  const std::vector<
    std::vector<RankingTree::TreeVertexIndex>
  >& treeRankingSets
) const {
  return temple::map(
    treeRankingSets,
    [&](const auto& set) -> std::vector<AtomIndexType> {
      return temple::map(
        set,
        [&](const auto& treeVertex) -> AtomIndexType {
          return _tree[treeVertex].molIndex;
        }
      );
    }
  );
}

std::vector<
  std::vector<AtomIndexType>
> RankingTree::getRanked() const {
  // We must transform the ranked tree vertex indices back to molecule indices.
  return _mapToAtomIndices(
    _branchOrderingHelper.getSets()
  );
}


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


} // namespace molassembler
