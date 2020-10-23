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
namespace Subgraphs {
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
    mappingsRef.get().push_back(makeIndexMap(a));

    //! Don't force stop after any mappings are found
    return true;
  }
};

template<typename Enum>
auto underlying(Enum a) {
  return static_cast<std::underlying_type_t<Enum>>(a);
}

struct PartialMolecule {
  using MaybePermutators = boost::optional<const StereopermutatorList&>;

  explicit PartialMolecule(const Molecule& molecule) : graph(molecule.graph()), permutatorsOption(molecule.stereopermutators()) {}
  explicit PartialMolecule(const Graph& graph) : graph(graph), permutatorsOption(boost::none) {}

  const Graph& graph;
  MaybePermutators permutatorsOption;
};

struct VertexComparator {
  VertexComparator(PartialMolecule a, PartialMolecule b, const VertexStrictness passStrictness)
    : data({std::move(a), std::move(b)}), strictness(passStrictness)
    {}

  bool operator () (const AtomIndex i, const AtomIndex j) const {
    if(data.first.graph.elementType(i) != data.second.graph.elementType(j)) {
      return false;
    }

    if(underlying(strictness) >= underlying(VertexStrictness::SubsumeShape)) {
      throw std::logic_error("Not implemented!");
      /* Determine if there is a low-effort (< 0.2 per step) transition from
       * the lower shape to the larger one, provided both vertices are
       * non-terminal in each molecule and both have a stereocenter defined on
       * it
       */
    }

    if(underlying(strictness) >= underlying(VertexStrictness::SubsumeStereopermutation)) {
      /* Determine if the shape on both vertices is the same, provided both
       * vertices are non-terminal
       */
      throw std::logic_error("Not implemented!");
    }

    return true;
  }

  std::pair<PartialMolecule, PartialMolecule> data;
  const VertexStrictness strictness;
};

/**
 * @brief Helper object to compare edges between two molecules
 */
struct EdgeComparator {
  EdgeComparator(PartialMolecule a, PartialMolecule b, const EdgeStrictness passStrictness)
    : data({std::move(a), std::move(b)}), strictness(passStrictness)
    {}

  bool operator () (const PrivateGraph::Edge i, const PrivateGraph::Edge j) const {
    if(
      underlying(strictness) >= underlying(EdgeStrictness::BondType)
      && data.first.graph.inner().bondType(i) != data.second.graph.inner().bondType(j)
    ) {
      return false;
    }

    if(underlying(strictness) >= underlying(EdgeStrictness::SubsumeStereopermutation)) {
      throw std::logic_error("Not implemented!");
    }

    return true;
  }

  std::pair<PartialMolecule, PartialMolecule> data;
  const EdgeStrictness strictness;
};

std::vector<IndexMap> maximumImpl(
  const PartialMolecule& a,
  const PartialMolecule& b,
  const VertexStrictness vertexStrictness,
  const EdgeStrictness edgeStrictness,
  const bool removeHydrogenPermutations
) {
  std::vector<IndexMap> mappings;
  SubgraphCallback callback {a.graph, mappings};

  boost::mcgregor_common_subgraphs_maximum_unique(
    a.graph.inner().bgl(),
    b.graph.inner().bgl(),
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
                  b.graph.elementType(fixedPair.second) == Utils::ElementType::H
                  && b.graph.elementType(sites.second) == Utils::ElementType::H
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

} // namespace

std::vector<IndexMap> maximum(
  const Graph& a,
  const Graph& b,
  const VertexStrictness vertexStrictness,
  const EdgeStrictness edgeStrictness,
  const bool removeHydrogenPermutations
) {
  if(underlying(vertexStrictness) >= underlying(VertexStrictness::SubsumeShape)) {
    throw std::runtime_error("Requested vertex comparison strictness not possible without stereopermutator information");
  }

  if(underlying(edgeStrictness) >= underlying(EdgeStrictness::SubsumeStereopermutation)) {
    throw std::runtime_error("Requested edge comparison strictness not possible without stereopermutator information");
  }

  return maximumImpl(
    PartialMolecule(a),
    PartialMolecule(b),
    vertexStrictness,
    edgeStrictness,
    removeHydrogenPermutations
  );
}

std::vector<IndexMap> maximum(
  const Molecule& a,
  const Molecule& b,
  const VertexStrictness vertexStrictness,
  const EdgeStrictness edgeStrictness,
  const bool removeHydrogenPermutations
) {
  return maximumImpl(
    PartialMolecule(a),
    PartialMolecule(b),
    vertexStrictness,
    edgeStrictness,
    removeHydrogenPermutations
  );
}

} // namespace Subgraphs
} // namespace Molassembler
} // namespace Scine
