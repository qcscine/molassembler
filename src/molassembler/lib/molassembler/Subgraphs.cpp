/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Subgraphs.h"

#include "boost/graph/mcgregor_common_subgraphs.hpp"

#include "shapes/PropertyCaching.h"

#include "shapes/Shapes.h"
#include "temple/constexpr/UpperTriangularMatrix.h"


#include "molassembler/Graph/Bridge.h"
#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"

// TEMP
// #include "Utils/Geometry/ElementInfo.h"
// #include <iostream>

/* TODO
 * - Missing algorithm for stereopermutator extension
 * - Ensure comparisons are symmetric! I.e. not needle-haystack but
 *   same-footing arguments
 */

namespace Scine {

namespace molassembler {

namespace subgraphs {

namespace detail {

struct SubgraphCallback {
  AtomIndex N;
  using IndexMapVector = std::vector<IndexMap>;
  std::reference_wrapper<IndexMapVector> mappingsRef;

  SubgraphCallback(
    const OuterGraph& a,
    IndexMapVector& mappings
  ) : N {a.N()},
      mappingsRef(mappings)
  {}

  template<class ABMap>
  IndexMap makeIndexMap(const ABMap& m) {
    IndexMap bimap;

    for(AtomIndex i = 0; i < N; ++i) {
      AtomIndex t = boost::get(m, i);
      if(t != boost::graph_traits<InnerGraph::BGLType>::null_vertex()) {
        bimap.insert(
          IndexMap::value_type(i, t)
        );
      }
    }

    return bimap;
  }

  template<class ABMap, class BAMap>
  bool operator() (
    ABMap a,
    BAMap /* b */,
    AtomIndex /* subgraphSize */
  ) {
    mappingsRef.get().push_back(
      makeIndexMap(a)
    );

    //! Don't force stop after any mappings are found
    return true;
  }
};

template<typename Enum>
auto underlying(Enum a) {
  return static_cast<
    std::underlying_type_t<Enum>
  >(a);
}

struct VertexComparator {
  using EnumUnderlying = std::underlying_type<VertexStrictness>::type;

  const Molecule& a;
  const Molecule& b;
  const VertexStrictness strictness;

  /* Maybe it is preferable to populate a cache of all combinations (w/ upper
   * triangular matrix, perhaps) the first time this is needed, and then just
   * use lookups
   */
  struct LowEffortMappingCache {
    static constexpr std::size_t upperTrigSize = Symmetry::nShapes * (Symmetry::nShapes - 1) / 2;

    static std::array<bool, upperTrigSize> defaultMatrixData() {
      std::array<bool, upperTrigSize> data;
      data.fill(false);
      return data;
    }

    //! Constructor establishing the matrix elements
    LowEffortMappingCache() : mappingMatrix(defaultMatrixData()) {
      // Populate with immediately adjacent information from symmetry mapping cache
      for(const auto from : Symmetry::allShapes) {
        for(const auto to : Symmetry::allShapes) {
          /* Smaller and equal sized to symmetries are always false and not stored
           * in the matrix at all. We can only anser directly adjacent cases,
           * so we skip everything else.
           */
          if(Symmetry::size(from) + 1 != Symmetry::size(to)) {
            continue;
          }

          // Get a reference to the matrix entry
          auto& matrixEntry = mappingMatrix.at(
            underlying(from),
            underlying(to)
          );

          auto mappingOptional = Symmetry::getMapping(from, to);

          // These should always be a Some type
          assert(mappingOptional);

          auto& mapping = mappingOptional.value();

          matrixEntry = (mapping.angularDistortion + mapping.chiralDistortion <= 0.2);
        }
      }

      // TODO!
      /* Now only directly adjacent transferability edges are present in the
       * matrix, but we want the matrix to contain more distant matching too
       */

      /* auto copyIn = [

      // Add transferability information
      for(const auto from : Symmetry::allShapes) {
        for(const auto to : Symmetry::allShapes) {
          // Skip smaller and equal-size target symmetries
          if(Symmetry::size(to) <= Symmetry::size(from)) {
            continue;
          }


        }
      }*/
    }

    bool subsumes(const Symmetry::Shape from, const Symmetry::Shape to) const {
      if(from == to) {
        return true;
      }

      return mappingMatrix.at(
        underlying(from),
        underlying(to)
      );
    }

    temple::UpperTriangularMatrix<bool, upperTrigSize> mappingMatrix;
  };

  static bool lowEffortMapping(
    const Symmetry::Shape from,
    const Symmetry::Shape to
  ) {
    static LowEffortMappingCache subsumptionMatrix;
    return subsumptionMatrix.subsumes(from, to);
  }

  VertexComparator(
    const Molecule& passA,
    const Molecule& passB,
    const VertexStrictness passStrictness
  ) : a(passA), b(passB), strictness(passStrictness) {}

  bool operator () (const AtomIndex i, const AtomIndex j) const {
    if(a.graph().elementType(i) != b.graph().elementType(j)) {
      return false;
    }

    if(underlying(strictness) >= underlying(VertexStrictness::SubsumeShape)) {
      throw std::logic_error("Not implemented!");
      /* Determine if there is a low-effort (< 0.2 per step) transition from
       * the lower shape to the larger one, provided both vertices are
       * non-terminal in each molecule and both have a stereocenter defined on
       * it
       */
      auto iAtomStereocenterOption = a.stereopermutators().option(i);
      auto jAtomStereocenterOption = b.stereopermutators().option(j);

      if(iAtomStereocenterOption && jAtomStereocenterOption) {
        Symmetry::Shape iShape = iAtomStereocenterOption->getShape();
        Symmetry::Shape jShape = jAtomStereocenterOption->getShape();

        if(iShape != jShape) {
          // Establish ordering for the call to lowEffortMapping
          if(Symmetry::size(iShape) > Symmetry::size(jShape)) {
            std::swap(iShape, jShape);
          }

          if(!lowEffortMapping(iShape, jShape)) {
            return false;
          }
        }
      }
    }

    if(underlying(strictness) >= underlying(VertexStrictness::SubsumeStereopermutation)) {
      /* Determine if the shape on both vertices is the same, provided both
       * vertices are non-terminal
       */
      throw std::logic_error("Not implemented!");
    }

    /*std::cout << "Matched " << Utils::ElementInfo::symbol(
      a.graph().elementType(i)
    ) << i << " to " << Utils::ElementInfo::symbol(
      b.graph().elementType(j)
    ) << j << "\n";*/

    return true;
  }
};

/**
 * @brief Helper object to compare edges between two molecules
 */
struct EdgeComparator {
  const Molecule& a;
  const Molecule& b;
  const EdgeStrictness strictness;

  EdgeComparator(
    const Molecule& passA,
    const Molecule& passB,
    EdgeStrictness passStrictness
  ) : a(passA), b(passB), strictness(passStrictness) {}

  bool operator () (const InnerGraph::Edge i, const InnerGraph::Edge j) const {
    /*auto bondStr = [](const InnerGraph& g, const InnerGraph::Edge& e) -> std::string {
      AtomIndex source = g.source(e);
      AtomIndex target = g.target(e);

      return (
        Utils::ElementInfo::symbol(g.elementType(source)) + std::to_string(source)
        + " - "
        + Utils::ElementInfo::symbol(g.elementType(target)) + std::to_string(target)
      );
    };
    std::cout << "Comparing " << bondStr(a.graph().inner(), i) << + " to " << bondStr(b.graph().inner(), j) << " -> ";*/

    if(
      underlying(strictness) >= underlying(EdgeStrictness::BondType)
      && a.graph().inner().bondType(i) != b.graph().inner().bondType(j)
    ) {
      return false;
    }

    if(underlying(strictness) >= underlying(EdgeStrictness::SubsumeStereopermutation)) {
      throw std::logic_error("Not implemented!");

      auto aBondStereocenterOptional = a.stereopermutators().option(
        toOuter(i, a.graph().inner())
      );
      auto bBondStereocenterOptional = b.stereopermutators().option(
        toOuter(j, b.graph().inner())
      );

      // Only if both have a stereocenter do we compare
      if(aBondStereocenterOptional && bBondStereocenterOptional) {
        // If both are assigned and have different values, the edges do not match
        if(
          aBondStereocenterOptional->assigned()
          && bBondStereocenterOptional->assigned()
          && aBondStereocenterOptional->numStereopermutations() == bBondStereocenterOptional->numStereopermutations()
          && (
            aBondStereocenterOptional->assigned().value()
            != bBondStereocenterOptional->assigned().value()
          )
        ) {
          return false;
        }
      }
    }

    return true;
  }
};

} // namespace detail

std::vector<IndexMap> maximum(
  const Molecule& a,
  const Molecule& b,
  const VertexStrictness vertexStrictness,
  const EdgeStrictness edgeStrictness,
  const bool removeHydrogenPermutations
) {
  std::vector<IndexMap> mappings;
  detail::SubgraphCallback callback {a.graph(), mappings};

  boost::mcgregor_common_subgraphs_maximum_unique(
    a.graph().inner().bgl(),
    b.graph().inner().bgl(),
    false, // Permit disconnected subgraphs
    callback,
    boost::vertices_equivalent(
      detail::VertexComparator {a, b, vertexStrictness}
    ).
    edges_equivalent(
      detail::EdgeComparator {a, b, edgeStrictness}
    )
  );

  if(removeHydrogenPermutations) {
    for(unsigned i = 0; i < mappings.size() - 1; ++i) {
      const IndexMap& fixedMap = mappings.at(i);
      mappings.erase(
        std::remove_if(
          std::begin(mappings) + i + 1,
          std::end(mappings),
          [&](const IndexMap& indexMap) -> bool {
            /* Left maps are ordered in their first index, so we can
             * sequentially compare elements of the mapping.
             *
             * The idea here is that if the mapping is sequence identical
             * regarding non-hydrogen elements, then it must be a permutation
             * thereof. This is much cheaper than using std::is_permutation.
             */
            return std::equal(
              std::begin(fixedMap.left),
              std::end(fixedMap.left),
              std::begin(indexMap.left),
              std::end(indexMap.left),
              [&](const auto& fixedPair, const auto& indexPair) -> bool {
                if(fixedPair.first != indexPair.first) {
                  return false;
                }

                if(
                  b.graph().elementType(fixedPair.second) == Utils::ElementType::H
                  && b.graph().elementType(indexPair.second) == Utils::ElementType::H
                ) {
                  return true;
                }

                return (fixedPair.second == indexPair.second);
              }
            );
          }
        ),
        std::end(mappings)
      );
    }
  }

  return mappings;
}

} // namespace subgraphs

} // namespace molassembler

} // namespace Scine
