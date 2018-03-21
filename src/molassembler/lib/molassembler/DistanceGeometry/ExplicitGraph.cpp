#include "DistanceGeometry/ExplicitGraph.h"

#include "Molecule.h"
#include "boost/graph/bellman_ford_shortest_paths.hpp"
#include "boost/graph/two_bit_color_map.hpp"
#include "DistanceGeometry/DistanceGeometry.h"
#include "DistanceGeometry/Error.h"

//#define USE_SPECIALIZED_GOR1_ALGORITHM
#ifdef USE_SPECIALIZED_GOR1_ALGORITHM
#include "DistanceGeometry/Gor1.h"
#else
#include "gor1/Gor1.h"
#endif

#include "temple/Random.h"

#include "AtomInfo.h"

/* Using Dijkstra's shortest paths despite there being negative edge weights is
 * alright since there are, by construction, no negative edge weight sum cycles,
 * which would be the death of the algorithm.
 *
 * But BGL does not allow this, so we would have to use Bellman-Ford, which is 
 * O(VE), not O(V + E).
 */

namespace molassembler {

namespace DistanceGeometry {

ExplicitGraph::ExplicitGraph(
  const Molecule& molecule,
  const BoundList& bounds
) : _graph {2 * molecule.numAtoms()},
    _molecule {molecule}
{
  // Populate the graph with bounds from list
  VertexDescriptor a, b;
  ValueBounds bound;
  for(const auto& boundTuple : bounds) {
    std::tie(a, b, bound) = boundTuple;

    addBound(a, b, bound);
  }

  addImplicitEdges();

  // Determine the two heaviest element types in the molecule, O(N)
  _heaviestAtoms = {{Delib::ElementType::H, Delib::ElementType::H}};
  unsigned N = molecule.numAtoms();
  for(AtomIndexType i = 0; i < N; ++i) {
    auto elementType = molecule.getElementType(i);
    if(
      static_cast<unsigned>(elementType) 
      > static_cast<unsigned>(_heaviestAtoms.back())
    ) {
      _heaviestAtoms.back() = elementType;

      if(
        static_cast<unsigned>(_heaviestAtoms.back()) 
        > static_cast<unsigned>(_heaviestAtoms.front())
      ) {
        std::swap(_heaviestAtoms.front(), _heaviestAtoms.back());
      }
    }
  }
}

void ExplicitGraph::_updateOrAddEdge(
  const VertexDescriptor& a,
  const VertexDescriptor& b,
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

  // The graph contains the lower bound negated
  return -boost::get(boost::edge_weight, _graph, edgeSearchPair.first);
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

double ExplicitGraph::maximalImplicitLowerBound(const VertexDescriptor& i) const {
  assert(isLeft(i));
  auto a = i % 2;
  auto elementType = _molecule.getElementType(a);

  if(elementType == _heaviestAtoms.front()) {
    return AtomInfo::vdwRadius(
      _heaviestAtoms.back()
    ) + AtomInfo::vdwRadius(elementType);
  }

  return AtomInfo::vdwRadius(
    _heaviestAtoms.front()
  ) + AtomInfo::vdwRadius(elementType);
}

const ExplicitGraph::GraphType& ExplicitGraph::getGraph() const {
  return _graph;
}

// This is O(N²)
outcome::result<Eigen::MatrixXd> ExplicitGraph::makeDistanceBounds() const noexcept {
  unsigned N = _molecule.numAtoms();

  Eigen::MatrixXd bounds;
  bounds.resize(N, N);
  bounds.setZero();

  unsigned M = boost::num_vertices(_graph);
  std::vector<double> distances (M);
  std::vector<VertexDescriptor> predecessors (M);
  using ColorMapType = boost::two_bit_color_map<>;
  ColorMapType color_map {M};

  for(AtomIndexType a = 0; a < N - 1; ++a) {
    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      boost::get(boost::vertex_index, _graph)
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      boost::get(boost::vertex_index, _graph)
    );

    // re-fill color map with white
    std::fill(
      color_map.data.get(),
      color_map.data.get() + (color_map.n + color_map.elements_per_char - 1) 
        / color_map.elements_per_char,
      0
    );

#ifdef USE_SPECIALIZED_GOR1_ALGORITHM
    boost::gor1_ig_shortest_paths(
      *this,
      VertexDescriptor {left(a)},
      predecessor_map,
      color_map,
      distance_map
    );
#else
    boost::gor1_simplified_shortest_paths(
      _graph,
      VertexDescriptor {left(a)},
      predecessor_map,
      color_map,
      distance_map
    );
#endif

    for(AtomIndexType b = a + 1; b < N; ++b) {
      bounds(a, b) = distances.at(left(b));
      bounds(b, a) = -distances.at(right(b));

      if(
        bounds(a, b) < bounds(b, a)
        || bounds(a, b) <= 0
        || bounds(b, a) <= 0
      ) {
        return DGError::GraphImpossible;
      }
    }
  }

  return bounds;
}

outcome::result<Eigen::MatrixXd> ExplicitGraph::makeDistanceMatrix() noexcept {
  return makeDistanceMatrix(Partiality::All);
}

outcome::result<Eigen::MatrixXd> ExplicitGraph::makeDistanceMatrix(Partiality partiality) noexcept {
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
    temple::random.randomEngine
  );

  unsigned M = boost::num_vertices(_graph);
  std::vector<double> distance (M);
  boost::two_bit_color_map<> color_map {M};
  std::vector<VertexDescriptor> predecessors (M);

  std::vector<AtomIndexType>::const_iterator separator;

  if(partiality == Partiality::FourAtom) {
    separator = indices.cbegin() + std::min(N, 4u);
  } else if(partiality == Partiality::TenPercent) {
    separator = indices.cbegin() + std::min(N, static_cast<unsigned>(0.1 * N));
  } else { // All
    separator = indices.cend();
  }

  for(auto iter = indices.cbegin(); iter != separator; ++iter) {
    const AtomIndexType a = *iter;

    std::vector<AtomIndexType> otherIndices;
    otherIndices.reserve(N - 1);

    // Avoid already-chosen elements
    for(AtomIndexType b = 0; b < a; ++b) {
      if(upperTriangle(b, a) == 0) {
        otherIndices.push_back(b);
      }
    }

    for(AtomIndexType b = a + 1; b < N; ++b) {
      if(upperTriangle(a, b) == 0) {
        otherIndices.push_back(b);
      }
    }

    std::shuffle(
      otherIndices.begin(),
      otherIndices.end(),
      temple::random.randomEngine
    );

    // Again through N - 1 indices: N²
    for(const auto& b : otherIndices) {

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

#ifdef USE_SPECIALIZED_GOR1_ALGORITHM
      boost::gor1_eg_shortest_paths(
        *this,
        left(a),
        predecessor_map,
        color_map,
        distance_map
      );

      double presumedLower = -distance.at(right(b));
      double presumedUpper = distance.at(left(b));

      double tightenedBound = temple::random.getSingle<double>(
        std::min(presumedLower, presumedUpper),
        std::max(presumedLower, presumedUpper)
      );
#else
      boost::gor1_simplified_shortest_paths(
        _graph,
        left(a),
        predecessor_map,
        color_map,
        distance_map
      );

      if(distance.at(left(b)) < -distance.at(right(b))) {
        return DGError::GraphImpossible;
      }

      /* Shortest distance from left a vertex to right b vertex is lower bound (negative)
       * Shortest distance from left a vertex to left b vertex is upper bound
       */
      double tightenedBound = temple::random.getSingle<double>(
        -distance.at(right(b)),
        distance.at(left(b))
      );
#endif

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

  for(auto iter = separator; iter != indices.cend(); ++iter) {
    const AtomIndexType a = *iter;

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

#ifdef USE_SPECIALIZED_GOR1_ALGORITHM
    boost::gor1_eg_shortest_paths(
      *this,
      left(a),
      predecessor_map,
      color_map,
      distance_map
    );
#else
    boost::gor1_simplified_shortest_paths(
      _graph,
      left(a),
      predecessor_map,
      color_map,
      distance_map
    );
#endif

    for(AtomIndexType b = 0; b < N; ++b) {
      if(a == b || upperTriangle(std::min(a, b), std::max(a, b)) > 0) {
        continue;
      }

      double presumedLower = -distance.at(right(b));
      double presumedUpper = distance.at(left(b));

      /* Shortest distance from left a vertex to right b vertex is lower bound (negative)
       * Shortest distance from left a vertex to left b vertex is upper bound
       */
      double tightenedBound = temple::random.getSingle<double>(
        std::min(presumedLower, presumedUpper),
        std::max(presumedLower, presumedUpper)
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

} // namespace molassembler
