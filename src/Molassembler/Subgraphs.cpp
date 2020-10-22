/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Subgraphs.h"

#include "boost/graph/mcgregor_common_subgraphs.hpp"

#include "Molassembler/Shapes/PropertyCaching.h"
#include "Molassembler/Shapes/Shapes.h"
#include "Molassembler/Temple/constexpr/UpperTriangularMatrix.h"
#include "Molassembler/Graph/Bridge.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/StereopermutatorList.h"

/* TODO
 * - Missing algorithm for stereopermutator extension
 * - Ensure comparisons are symmetric! I.e. not needle-haystack but
 *   same-footing arguments
 */

namespace Scine {
namespace Molassembler {
namespace subgraphs {
namespace {

struct SubgraphCallback {
  AtomIndex N;
  using IndexMapVector = std::vector<IndexMap>;
  std::reference_wrapper<IndexMapVector> mappingsRef;

  SubgraphCallback(
    const Graph& a,
    IndexMapVector& mappings
  ) : N {a.N()},
      mappingsRef(mappings)
  {}

  template<class AbMap>
  IndexMap makeIndexMap(const AbMap& m) {
    IndexMap bimap;

    for(AtomIndex i = 0; i < N; ++i) {
      AtomIndex t = boost::get(m, i);
      if(t != boost::graph_traits<PrivateGraph::BglType>::null_vertex()) {
        bimap.insert(
          IndexMap::value_type(i, t)
        );
      }
    }

    return bimap;
  }

  template<class AbMap, class BaMap>
  bool operator() (
    AbMap a,
    BaMap /* b */,
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
    static constexpr std::size_t upperTrigSize = Shapes::nShapes * (Shapes::nShapes - 1) / 2;

    static std::array<bool, upperTrigSize> defaultMatrixData() {
      std::array<bool, upperTrigSize> data;
      data.fill(false);
      return data;
    }

    //! Constructor establishing the matrix elements
    LowEffortMappingCache() : mappingMatrix(defaultMatrixData()) {
      // Populate with immediately adjacent information from symmetry mapping cache
      for(const auto from : Shapes::allShapes) {
        for(const auto to : Shapes::allShapes) {
          /* Smaller and equal sized to symmetries are always false and not stored
           * in the matrix at all. We can only anser directly adjacent cases,
           * so we skip everything else.
           */
          if(Shapes::size(from) + 1 != Shapes::size(to)) {
            continue;
          }

          // Get a reference to the matrix entry
          auto& matrixEntry = mappingMatrix.at(
            underlying(from),
            underlying(to)
          );

          auto mappingOptional = Shapes::getMapping(from, to);

          // These should always be a Some type
          assert(mappingOptional);

          const auto& mapping = mappingOptional.value();

          matrixEntry = (mapping.angularDistortion + mapping.chiralDistortion <= 0.2);
        }
      }

      // TODO!
      /* Now only directly adjacent transferability edges are present in the
       * matrix, but we want the matrix to contain more distant matching too
       */

      /* auto copyIn = [

      // Add transferability information
      for(const auto from : Shapes::allShapes) {
        for(const auto to : Shapes::allShapes) {
          // Skip smaller and equal-size target symmetries
          if(Shapes::size(to) <= Shapes::size(from)) {
            continue;
          }


        }
      }*/
    }

    bool subsumes(const Shapes::Shape from, const Shapes::Shape to) const {
      if(from == to) {
        return true;
      }

      return mappingMatrix.at(
        underlying(from),
        underlying(to)
      );
    }

    Temple::UpperTriangularMatrix<bool, upperTrigSize> mappingMatrix;
  };

  static bool lowEffortMapping(
    const Shapes::Shape from,
    const Shapes::Shape to
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
        Shapes::Shape iShape = iAtomStereocenterOption->getShape();
        Shapes::Shape jShape = jAtomStereocenterOption->getShape();

        if(iShape != jShape) {
          // Establish ordering for the call to lowEffortMapping
          if(Shapes::size(iShape) > Shapes::size(jShape)) {
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

  bool operator () (const PrivateGraph::Edge i, const PrivateGraph::Edge j) const {
    /*auto bondStr = [](const PrivateGraph& g, const PrivateGraph::Edge& e) -> std::string {
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

} // namespace

std::vector<IndexMap> maximum(
  const Molecule& a,
  const Molecule& b,
  const VertexStrictness vertexStrictness,
  const EdgeStrictness edgeStrictness,
  const bool removeHydrogenPermutations
) {
  std::vector<IndexMap> mappings;
  SubgraphCallback callback {a.graph(), mappings};

  boost::mcgregor_common_subgraphs_maximum_unique(
    a.graph().inner().bgl(),
    b.graph().inner().bgl(),
    false, // Permit disconnected subgraphs
    callback,
    boost::vertices_equivalent(
      VertexComparator {a, b, vertexStrictness}
    ).
    edges_equivalent(
      EdgeComparator {a, b, edgeStrictness}
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
              [&](const auto& fixedPair, const auto& sites) -> bool {
                if(fixedPair.first != sites.first) {
                  return false;
                }

                if(
                  b.graph().elementType(fixedPair.second) == Utils::ElementType::H
                  && b.graph().elementType(sites.second) == Utils::ElementType::H
                ) {
                  return true;
                }

                return (fixedPair.second == sites.second);
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
} // namespace Molassembler
} // namespace Scine
