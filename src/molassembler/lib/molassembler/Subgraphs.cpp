#include "molassembler/Subgraphs.h"

#include "boost/graph/mcgregor_common_subgraphs.hpp"

#include "chemical_symmetries/Properties.h"

#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/StereocenterList.h"

namespace molassembler {

namespace subgraphs {

namespace detail {

struct SubgraphCallback {
  AtomIndex N;
  std::vector<IndexMap> mappings;

  SubgraphCallback() = default;
  SubgraphCallback(const OuterGraph& a) : N {a.N()} {}

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
    mappings.push_back(
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

  /* Maybe its preferable to populate a cache of all combinations (w/ upper
   * triangular matrix, perhaps) the first time this is needed, and then just
   * use lookups
   */

  static bool lowEffortMapping(
    const Symmetry::Name small,
    const Symmetry::Name large
  ) {
    const unsigned largeSize = Symmetry::size(large);
    unsigned currentSize = Symmetry::size(small);

    std::vector<Symmetry::Name> currentSymmetries {small};
    while(largeSize > currentSize + 1 && !currentSymmetries.empty()) {
      std::vector<Symmetry::Name> nextSizeSymmetries;

      /* Create a list of next-size symmetries that have low-effort mappings
       * from the current symmetries. Need to avoid duplicates, so perhaps use
       * TinySets? Do those have the required behavior of not reinserting if the
       * value is already present?
       */
    }

    if(currentSymmetries.empty()) {
      return false;
    }

    for(const auto symmetry : currentSymmetries) {
      /* If any of the current symmetries have a low-effort mapping to the large
       * symmetry, return true
       */
      // TODO continue
    }

    return false;
  }

  VertexComparator(
    const Molecule& passA,
    const Molecule& passB,
    VertexStrictness passStrictness
  ) : a(passA), b(passB), strictness(passStrictness) {}

  bool operator () (const AtomIndex i, const AtomIndex j) const {
    if(a.graph().elementType(i) != b.graph().elementType(j)) {
      return false;
    }

    if(
      underlying(strictness)
      >= underlying(VertexStrictness::LowEffortTransitionToLargerSymmetry)
    ) {
      /* Determine if there is a low-effort (< 0.2 per step) transition from
       * the lower symmetry to the larger one, provided both vertices are
       * non-terminal in each molecule and both have a stereocenter defined on
       * it
       */
      auto iAtomStereocenterOption = a.stereocenters().option(i);
      auto jAtomStereocenterOption = b.stereocenters().option(j);

      if(iAtomStereocenterOption && jAtomStereocenterOption) {
        Symmetry::Name iSymmetry = iAtomStereocenterOption->getSymmetry();
        Symmetry::Name jSymmetry = jAtomStereocenterOption->getSymmetry();

        if(iSymmetry != jSymmetry) {
          // Establish ordering for the call to lowEffortMapping
          if(Symmetry::size(iSymmetry) > Symmetry::size(jSymmetry)) {
            std::swap(iSymmetry, jSymmetry);
          }

          if(!lowEffortMapping(iSymmetry, jSymmetry)) {
            return false;
          }
        }
      }
    }

    if(underlying(strictness) >= underlying(VertexStrictness::SameSymmetry)) {
      /* Determine if the symmetry on both vertices is the same, provided both
       * vertices are non-terminal
       */
    }

    return true;
  }
};

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
    if(a.graph().inner().bondType(i) != b.graph().inner().bondType(j)) {
      return false;
    }

    if(strictness == EdgeStrictness::EZIdentical) {
      auto aBondStereocenterOptional = a.stereocenters().option(
        BondIndex {
          a.graph().inner().source(i),
          a.graph().inner().target(i)
        }
      );
      auto bBondStereocenterOptional = b.stereocenters().option(
        BondIndex {
          b.graph().inner().source(j),
          b.graph().inner().target(j)
        }
      );

      // Only if both have a stereocenter do we compare
      if(aBondStereocenterOptional && bBondStereocenterOptional) {
        // If both are assigned and have different values, the edges do not match
        if(
          aBondStereocenterOptional->assigned()
          && bBondStereocenterOptional->assigned()
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

std::unique_ptr<
  temple::UpperTriangularMatrix<bool, upperTrigSize>
> lowEffortMappings;

void generateLowEffortMappings() {
  lowEffortMappings = std::make_unique<
    temple::UpperTriangularMatrix<bool, upperTrigSize>
  >();

  // Populate with immediately adjacent information from symmetry mapping cache
  for(const auto from : Symmetry::allNames) {
    for(const auto to : Symmetry::allNames) {
      auto& matrixEntry = lowEffortMappings -> at(
        underlying(from),
        underlying(to)
      );

      // Smaller and equal sized to symmetries are false
      if(Symmetry::size(to) <= Symmetry::size(from)) {
        matrixEntry = false;
        continue;
      }

      // Initially, transitions to larger size than size(from) + 1 are false
      if(Symmetry::size(to) > Symmetry::size(from) + 1) {
        matrixEntry = false;
        continue;
      }

      // Remaining case: size(to) == from + 1
      auto mappingOptional = Symmetry::getMapping(from, to);

      // These should always be a Some type
      assert(mappingOptional);

      auto& mapping = mappingOptional.value();

      matrixEntry = (mapping.angularDistortion + mapping.chiralDistortion <= 0.2);
    }
  }

  /* auto copyIn = [

  // Add transferability information
  for(const auto from : Symmetry::allNames) {
    for(const auto to : Symmetry::allNames) {
      // Skip smaller and equal-size target symmetries
      if(Symmetry::size(to) <= Symmetry::size(from)) {
        continue;
      }


    }
  }*/
}

} // namespace detail

std::vector<IndexMap> maximum(
  const Molecule& a,
  const Molecule& b,
  VertexStrictness vertexStrictness,
  EdgeStrictness edgeStrictness
) {
  if(vertexStrictness == VertexStrictness::LowEffortTransitionToLargerSymmetry && !(detail::lowEffortMappings)) {
    detail::generateLowEffortMappings();
  }

  detail::SubgraphCallback callback;

  boost::mcgregor_common_subgraphs_maximum(
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

  return callback.mappings;
}

} // namespace subgraphs

} // namespace molassembler
