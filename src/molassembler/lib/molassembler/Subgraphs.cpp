#include "Subgraphs.h"

#include "boost/graph/mcgregor_common_subgraphs.hpp"

#include "GraphHelpers.h"

namespace molassembler {

namespace subgraphs {

struct SubgraphCallback {
  AtomIndexType N;
  std::vector<IndexMap> mappings;

  SubgraphCallback() = default;
  SubgraphCallback(const GraphType& a) : N {graph::numVertices(a)} {}

  template<class NeedleToHaystackMap>
  IndexMap makeIndexMap(const NeedleToHaystackMap& m) {
    IndexMap mapping;

    for(AtomIndexType i = 0; i < N; ++i) {
      AtomIndexType t = boost::get(m, i);
      if(t != boost::graph_traits<GraphType>::null_vertex()) {
        mapping.emplace(i, t);
      }
    }

    return mapping;
  }

  template<class NeedleToHaystackMap, class HaystackToNeedleMap>
  bool operator() (
    NeedleToHaystackMap a,
    HaystackToNeedleMap /* b */,
    AtomIndexType /* subgraphSize */
  ) {
    mappings.push_back(
      makeIndexMap(a)
    );

    //! Don't force stop after any mappings are found
    return true;
  }
};

std::vector<IndexMap> maximum(
  const GraphType& a,
  const GraphType& b
) {
  SubgraphCallback callback;

  boost::mcgregor_common_subgraphs_maximum(
    a,
    b,
    false, // Permit disconnected subgraphs
    callback,
    boost::vertices_equivalent(
      [&a, &b](const AtomIndexType i, const AtomIndexType j) -> bool {
        return a[i].elementType == b[j].elementType;
      }
    ).
    edges_equivalent(
      [&a, &b](const GraphType::edge_descriptor i, const GraphType::edge_descriptor j) -> bool {
        return a[i].bondType == b[j].bondType;
      }
    )
  );

  return callback.mappings;
}

std::vector<IndexMap> maximumUnique(
  const GraphType& a,
  const GraphType& b
) {
  SubgraphCallback callback;

  boost::mcgregor_common_subgraphs_maximum_unique(
    a,
    b,
    false, // Permit disconnected subgraphs
    callback,
    boost::vertices_equivalent(
      [&a, &b](const AtomIndexType i, const AtomIndexType j) -> bool {
        return a[i].elementType == b[j].elementType;
      }
    ).
    edges_equivalent(
      [&a, &b](const GraphType::edge_descriptor i, const GraphType::edge_descriptor j) -> bool {
        return a[i].bondType == b[j].bondType;
      }
    )
  );

  return callback.mappings;
}

} // namespace subgraphs

} // namespace molassembler
