#include "DistanceGeometry/ExplicitGraph.h"

#include "boost/graph/bellman_ford_shortest_paths.hpp"
#include "boost/graph/two_bit_color_map.hpp"
#include "Graph/Gor1.h"
#include "template_magic/Random.h"

#include "AtomInfo.h"

/* Using Dijkstra's shortest paths despite there being negative edge weights is
 * alright since there are, by construction, no negative edge weight sum cycles,
 * which would be the death of the algorithm.
 *
 * But BGL does not allow this, so we have to use Bellman-Ford, which is 
 * O(VE), not O(V + E).
 */

namespace MoleculeManip {

namespace DistanceGeometry {

ExplicitGraph::ExplicitGraph(
  const Molecule& molecule,
  const BoundList& bounds
) : _graph {2 * molecule.numAtoms()},
    _generatedDistances {false},
    _molecule {molecule}
{
  // Populate the graph with bounds from list
  ExplicitGraphType::vertex_descriptor a, b;
  ValueBounds bound;
  for(const auto& boundTuple : bounds) {
    std::tie(a, b, bound) = boundTuple;

    addBound(a, b, bound);
  }

  addImplicitEdges();
}

void ExplicitGraph::_updateOrAddEdge(
  const ExplicitGraphType::vertex_descriptor& a,
  const ExplicitGraphType::vertex_descriptor& b,
  const double& edgeWeight
) {
  auto edgeSearchPair = boost::edge(a, b, _graph);
  if(edgeSearchPair.second) {
    boost::get(boost::edge_weight, _graph, edgeSearchPair.first) = edgeWeight;
  } else {
    boost::add_edge(a, b, edgeWeight, _graph);
  }
}

void ExplicitGraph::_updateGraphWithFixedDistance(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const double& fixedDistance
) {
  _updateOrAddEdge(left(a), left(b), fixedDistance);
  _updateOrAddEdge(left(b), left(a), fixedDistance);

  _updateOrAddEdge(right(a), right(b), fixedDistance);
  _updateOrAddEdge(right(b), right(a), fixedDistance);

  _updateOrAddEdge(left(a), right(b), -fixedDistance);
  _updateOrAddEdge(left(b), right(a), -fixedDistance);
}

void ExplicitGraph::addBound(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const ValueBounds& bound
) {
  // Bidirectional edge in left graph with upper weight
  boost::add_edge(left(a), left(b), bound.upper, _graph);
  boost::add_edge(left(b), left(a), bound.upper, _graph);

  // Bidirectional edge in right graph with upper weight
  boost::add_edge(right(a), right(b), bound.upper, _graph);
  boost::add_edge(right(b), right(a), bound.upper, _graph);

  // Forward edge from left to right graph with negative lower bound weight
  boost::add_edge(left(a), right(b), -bound.lower, _graph);
  boost::add_edge(left(b), right(a), -bound.lower, _graph);
}

void ExplicitGraph::addImplicitEdges() {
  AtomIndexType N = boost::num_vertices(_graph) / 2;

  for(AtomIndexType a = 0; a < N; ++a) {
    for(AtomIndexType b = a + 1; b < N; ++b) {
      double vdwLowerBound = -1 * (
        AtomInfo::vdwRadius(
          _molecule.getElementType(a)
        ) + AtomInfo::vdwRadius(
          _molecule.getElementType(b)
        )
      );

      auto edgeSearchPair = boost::edge(left(a), right(b), _graph);
      if(!edgeSearchPair.second) {
        boost::add_edge(left(a), right(b), vdwLowerBound, _graph);
      }

      edgeSearchPair = boost::edge(left(b), right(a), _graph);
      if(!edgeSearchPair.second) {
        boost::add_edge(left(b), right(a), vdwLowerBound, _graph);
      }
    }
  }
}

double ExplicitGraph::lowerBound(
  const AtomIndexType& i,
  const AtomIndexType& j
) const {
  auto edgeSearchPair = boost::edge(
    left(i),
    right(j),
    _graph
  );

  assert(edgeSearchPair.second);

  return boost::get(boost::edge_weight, _graph, edgeSearchPair.first);
}

double ExplicitGraph::upperBound(
  const AtomIndexType& i,
  const AtomIndexType& j
) const {
  auto edgeSearchPair = boost::edge(
    left(i),
    left(j),
    _graph
  );

  assert(edgeSearchPair.second);

  return boost::get(boost::edge_weight, _graph, edgeSearchPair.first);
}

const ExplicitGraph::ExplicitGraphType& ExplicitGraph::getGraph() const {
  return _graph;
}

// This is O(N²)
DistanceBoundsMatrix ExplicitGraph::makeDistanceBounds() const {
  unsigned N = _molecule.numAtoms();

  DistanceBoundsMatrix bounds {N};

  for(AtomIndexType i = 0; i < N - 1; ++i) {
    unsigned M = boost::num_vertices(_graph);

    std::vector<double> distance (M);

    boost::bellman_ford_shortest_paths(
      _graph,
      M,
      boost::root_vertex(left(i)).
      distance_map(&distance[0])
    );

    for(AtomIndexType j = i + 1; j < N; ++j) {
      /* Since the edge weights from left to right are negative, the lower
       * bounds are negative too. Include the minimum lower bounds in case
       * the shortest path is greater (i.e. no explicit shortest path is
       * considered) than zero.
       */
      double lowerBound;

      if(distance.at(right(j)) >= 0) {
        lowerBound = AtomInfo::vdwRadius(
          _molecule.getElementType(i)
        ) + AtomInfo::vdwRadius(
          _molecule.getElementType(j)
        );
      } else {
        lowerBound = -distance.at(right(j));
      }

      assert(lowerBound > 0);
      assert(lowerBound < distance.at(left(j)));

      bounds.setLowerBound(i, j, lowerBound);
      bounds.setUpperBound(i, j, distance.at(left(j)));
    }
  }

  return bounds;
}

Eigen::MatrixXd ExplicitGraph::makeDistanceMatrix() {
  if(_generatedDistances) {
    throw std::logic_error("ExplicitGraph: Already generated distances destructively!");
  }

  _generatedDistances = true;

  unsigned N = _molecule.numAtoms();

  Eigen::MatrixXd distancesMatrix;
  distancesMatrix.resize(N, N);
  distancesMatrix.setZero();

  auto upperTriangle = distancesMatrix.triangularView<Eigen::StrictlyUpper>();

  std::vector<AtomIndexType> indices (N);

  std::iota(
    indices.begin(),
    indices.end(),
    0
  );

  std::shuffle(
    indices.begin(),
    indices.end(),
    TemplateMagic::random.randomEngine
  );

  unsigned M = boost::num_vertices(_graph);
  std::vector<double> distance (M);
  boost::two_bit_color_map<> color_map {M};
  std::vector<ExplicitGraphType::vertex_descriptor> predecessors (M);

  // Once through N indices: N
  for(const auto& a : indices) {
    std::vector<AtomIndexType> otherIndices (N - 1);
    std::iota(
      otherIndices.begin(),
      otherIndices.begin() + a,
      0
    );

    std::iota(
      otherIndices.begin() + a,
      otherIndices.end(),
      a
    );

    std::shuffle(
      otherIndices.begin(),
      otherIndices.end(),
      TemplateMagic::random.randomEngine
    );

    // Again through N - 1 indices: N²
    for(const auto& b : otherIndices) {
      if(
        a == b
        || upperTriangle(
          std::min(a, b),
          std::max(a, b)
        ) > 0
      ) {
        // Skip on-diagonal and already-chosen entries
        continue;
      }

      auto predecessor_map = boost::make_iterator_property_map(
        predecessors.begin(),
        boost::get(boost::vertex_index, _graph)
      );

      auto distance_map = boost::make_iterator_property_map(
        distance.begin(),
        boost::get(boost::vertex_index, _graph)
      );

      // re-fill color map with white
      std::fill(
        color_map.data.get(),
        color_map.data.get() + (color_map.n + color_map.elements_per_char - 1) 
          / color_map.elements_per_char,
        0
      );

      boost::gor1_simplified_shortest_paths(
        _graph,
        left(a),
        predecessor_map,
        color_map,
        distance_map
      );

      /* Since the edge weights from left to right are negative, the lower
       * bounds are negative too. Include the minimum lower bounds in case
       * the shortest path is greater (i.e. no explicit shortest path is
       * considered) than zero.
       */
      double lowerBound = -distance.at(right(b));

      assert(lowerBound > 0);
      assert(lowerBound < distance.at(left(b)));

      /* Shortest distance from left a vertex to right b vertex is lower bound (negative)
       * Shortest distance from left a vertex to left b vertex is upper bound
       */
      double tightenedBound = TemplateMagic::random.getSingle<double>(
        lowerBound,
        distance.at(left(b))
      );

      upperTriangle(
        std::min(a, b),
        std::max(a, b)
      ) = tightenedBound;

      // Modify the graph accordingly
      _updateGraphWithFixedDistance(
        a,
        b,
        tightenedBound
      );
    }
  }

  return distancesMatrix;
}

} // namespace DistanceGeometry

} // namespace MoleculeManip
