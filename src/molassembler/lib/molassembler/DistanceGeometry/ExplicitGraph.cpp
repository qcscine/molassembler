/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/DistanceGeometry/ExplicitGraph.h"

#include "boost/graph/two_bit_color_map.hpp"

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/DistanceGeometry/Error.h"
#include "molassembler/Log.h"
#include "molassembler/Modeling/AtomInfo.h"
#include "molassembler/Molecule.h"
#include "molassembler/Options.h"
#include "molassembler/Graph/InnerGraph.h"

#include "temple/Random.h"

#ifdef MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
#include "molassembler/DistanceGeometry/Gor1.h"
#else
#include "gor1/Gor1.h"
#endif


/* Using Dijkstra's shortest paths despite there being negative edge weights is
 * alright since there are, by construction, no negative edge weight sum cycles,
 * which would be the death of the algorithm.
 *
 * But BGL does not allow this, so we would have to use Bellman-Ford, which is
 * O(VE), not O(V + E). Instead of accepting this, we use GOR1 instead, which
 * is faster on graphs with many negative weights.
 */

namespace Scine {

namespace molassembler {

namespace DistanceGeometry {

ExplicitGraph::ExplicitGraph(
  const InnerGraph& inner,
  const BoundsMatrix& bounds
) : _graph {2 * inner.N()},
    _inner {inner}
{
  const AtomIndex N = inner.N();

  for(AtomIndex a = 0; a < N; ++a) {
    for(AtomIndex b = a + 1; b < N; ++b) {
      // a < b in all iterations
      const double lowerBound = bounds(b, a);
      const double upperBound = bounds(a, b);
      assert(lowerBound <= upperBound);

      if(lowerBound == 0.0 && upperBound == 0.0) {
        /* Add implicit edges for this a-b pair. "Implicit" edge (minimum
         * distance is sum of vdw radii) is explicit in this graph variant.
        */
        double vdwLowerBound = (
          AtomInfo::vdwRadius(inner.elementType(a))
          + AtomInfo::vdwRadius(inner.elementType(b))
        );

        boost::add_edge(left(a), right(b), -vdwLowerBound, _graph);
        boost::add_edge(left(b), right(a), -vdwLowerBound, _graph);
      } else {
        /* Add explicit edges for this a-b pair */
        // Bidirectional edge in left graph with upper weight
        boost::add_edge(left(a), left(b), upperBound, _graph);
        boost::add_edge(left(b), left(a), upperBound, _graph);

        // Bidirectional edge in right graph with upper weight
        boost::add_edge(right(a), right(b), upperBound, _graph);
        boost::add_edge(right(b), right(a), upperBound, _graph);

        // Forward edge from left to right graph with negative lower bound weight
        boost::add_edge(left(a), right(b), -lowerBound, _graph);
        boost::add_edge(left(b), right(a), -lowerBound, _graph);
      }
    }
  }

  // Determine the two heaviest element types in the molecule, O(N)
  _heaviestAtoms = {{Utils::ElementType::H, Utils::ElementType::H}};
  for(AtomIndex i = 0; i < N; ++i) {
    auto elementType = inner.elementType(i);
    if(
      Utils::ElementInfo::Z(elementType)
      > Utils::ElementInfo::Z(_heaviestAtoms.back())
    ) {
      _heaviestAtoms.back() = elementType;

      if(
        Utils::ElementInfo::Z(_heaviestAtoms.back())
        > Utils::ElementInfo::Z(_heaviestAtoms.front())
      ) {
        std::swap(_heaviestAtoms.front(), _heaviestAtoms.back());
      }
    }
  }
}

ExplicitGraph::ExplicitGraph(
  const InnerGraph& inner,
  const DistanceBoundsMatrix& bounds
) : _graph {2 * inner.N()},
    _inner {inner}
{
  const VertexDescriptor N = inner.N();
  for(VertexDescriptor a = 0; a < N; ++a) {
    for(VertexDescriptor b = a + 1; b < N; ++b) {
      const double lower = bounds.lowerBound(a, b);
      const double upper = bounds.upperBound(a, b);

      if(lower != DistanceBoundsMatrix::defaultLower) {
        // Forward edge from left to right graph with negative lower bound weight
        boost::add_edge(left(a), right(b), -lower, _graph);
        boost::add_edge(left(b), right(a), -lower, _graph);
      } else {
        const double vdwLowerBound = (
          AtomInfo::vdwRadius(_inner.elementType(a))
          + AtomInfo::vdwRadius(_inner.elementType(b))
        );

        // Implicit lower bound on distance between the vertices
        boost::add_edge(left(a), right(b), -vdwLowerBound, _graph);
        boost::add_edge(left(b), right(a), -vdwLowerBound, _graph);
      }

      if(upper != DistanceBoundsMatrix::defaultUpper) {
        // Bidirectional edge in left graph with upper weight
        boost::add_edge(left(a), left(b), upper, _graph);
        boost::add_edge(left(b), left(a), upper, _graph);

        // Bidirectional edge in right graph with upper weight
        boost::add_edge(right(a), right(b), upper, _graph);
        boost::add_edge(right(b), right(a), upper, _graph);
      }
    }
  }

  // Determine the two heaviest element types in the molecule, O(N)
  _heaviestAtoms = {{Utils::ElementType::H, Utils::ElementType::H}};
  for(AtomIndex i = 0; i < N; ++i) {
    auto elementType = inner.elementType(i);
    if(
      Utils::ElementInfo::Z(elementType)
      > Utils::ElementInfo::Z(_heaviestAtoms.back())
    ) {
      _heaviestAtoms.back() = elementType;

      if(
        Utils::ElementInfo::Z(_heaviestAtoms.back())
        > Utils::ElementInfo::Z(_heaviestAtoms.front())
      ) {
        std::swap(_heaviestAtoms.front(), _heaviestAtoms.back());
      }
    }
  }
}

void ExplicitGraph::addBound(
  const VertexDescriptor a,
  const VertexDescriptor b,
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

void ExplicitGraph::explainContradictionPaths(
  const VertexDescriptor a,
  const VertexDescriptor b,
  const std::vector<VertexDescriptor>& predecessors,
  const std::vector<double>& distances
) {
  using LevelBaseType = std::underlying_type<Log::Level>::type;
  if(
    static_cast<LevelBaseType>(Log::level)
    <= static_cast<LevelBaseType>(Log::Level::Warning)
  ) {
    // Report the contradiction in the log
    auto& logRef = Log::log(Log::Level::Warning);
    logRef << "Encountered contradiction in triangle inequality limits calculation.\n";
    logRef << "Path in graph for upper bound: l" << b;

    AtomIndex intermediate = left(b);
    do {
      intermediate = predecessors[intermediate];
      logRef << " <- l" << (intermediate / 2);
    } while (intermediate != left(a));
    logRef << ". Length " << distances.at(left(b))
      << "\nPath in graph for lower bound: r" << b;

    intermediate = right(b);
    do {
      intermediate = predecessors[intermediate];
      logRef << " <- "
        << (intermediate % 2 == 0 ? "l" : "r")
        << (intermediate / 2);
    } while (intermediate != left(a));

    logRef << ". Length " << distances.at(right(b)) << "\n";
  }
}

void ExplicitGraph::_updateOrAddEdge(
  const VertexDescriptor i,
  const VertexDescriptor j,
  const double edgeWeight
) {
  auto edgeSearchPair = boost::edge(i, j, _graph);
  if(edgeSearchPair.second) {
    boost::get(boost::edge_weight, _graph, edgeSearchPair.first) = edgeWeight;
  } else {
    boost::add_edge(i, j, edgeWeight, _graph);
  }
}

void ExplicitGraph::_updateGraphWithFixedDistance(
  const VertexDescriptor a,
  const VertexDescriptor b,
  const double fixedDistance
) {
  _updateOrAddEdge(left(a), left(b), fixedDistance);
  _updateOrAddEdge(left(b), left(a), fixedDistance);

  _updateOrAddEdge(right(a), right(b), fixedDistance);
  _updateOrAddEdge(right(b), right(a), fixedDistance);

  _updateOrAddEdge(left(a), right(b), -fixedDistance);
  _updateOrAddEdge(left(b), right(a), -fixedDistance);
}

double ExplicitGraph::lowerBound(
  const VertexDescriptor a,
  const VertexDescriptor b
) const {
  auto edgeSearchPair = boost::edge(left(a), right(b), _graph);

  assert(edgeSearchPair.second);

  // The graph contains the lower bound negated
  return -boost::get(boost::edge_weight, _graph, edgeSearchPair.first);
}

double ExplicitGraph::upperBound(
  const VertexDescriptor a,
  const VertexDescriptor b
) const {
  auto edgeSearchPair = boost::edge(left(a), left(b), _graph);

  assert(edgeSearchPair.second);

  return boost::get(boost::edge_weight, _graph, edgeSearchPair.first);
}

double ExplicitGraph::maximalImplicitLowerBound(const VertexDescriptor i) const {
  assert(isLeft(i));
  AtomIndex a = i / 2;
  Utils::ElementType elementType = _inner.elementType(a);

  if(elementType == _heaviestAtoms.front()) {
    return AtomInfo::vdwRadius(
      _heaviestAtoms.back()
    ) + AtomInfo::vdwRadius(elementType);
  }

  return AtomInfo::vdwRadius(
    _heaviestAtoms.front()
  ) + AtomInfo::vdwRadius(elementType);
}

const ExplicitGraph::GraphType& ExplicitGraph::graph() const {
  return _graph;
}

outcome::result<Eigen::MatrixXd> ExplicitGraph::makeDistanceBounds() const noexcept {
  unsigned N = _inner.N();

  Eigen::MatrixXd bounds;
  bounds.resize(N, N);
  bounds.setZero();

  unsigned M = boost::num_vertices(_graph);
  std::vector<double> distances (M);
  std::vector<VertexDescriptor> predecessors (M);
  using ColorMapType = boost::two_bit_color_map<>;
  ColorMapType color_map {M};

  for(AtomIndex a = 0; a < N - 1; ++a) {
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
      color_map.data.get() + (color_map.n + ColorMapType::elements_per_char - 1)
        / ColorMapType::elements_per_char,
      0
    );

#ifdef MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
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

    for(AtomIndex b = a + 1; b < N; ++b) {
      // Get upper bound from distances
      bounds(a, b) = distances.at(left(b));
      // Get lower bound from distances
      bounds(b, a) = -distances.at(right(b));

      // If the upper bound is less than the lower bound, we have a contradiction
      if(bounds(a, b) < bounds(b, a)) {
        explainContradictionPaths(a, b, predecessors, distances);
        return DgError::GraphImpossible;
      }

      // Negative values are not allowed
      if(bounds(a, b) <= 0 || bounds(b, a) <= 0) {
        return DgError::GraphImpossible;
      }
    }
  }

  return bounds;
}

outcome::result<Eigen::MatrixXd> ExplicitGraph::makeDistanceMatrix(random::Engine& engine) noexcept {
  return makeDistanceMatrix(engine, Partiality::All);
}

outcome::result<Eigen::MatrixXd> ExplicitGraph::makeDistanceMatrix(random::Engine& engine, Partiality partiality) noexcept {
  const unsigned N = _inner.N();

  Eigen::MatrixXd distancesMatrix;
  distancesMatrix.resize(N, N);
  distancesMatrix.setZero();

  auto upperTriangle = distancesMatrix.triangularView<Eigen::StrictlyUpper>();

  std::vector<AtomIndex> indices (N);

  std::iota(std::begin(indices), std::end(indices), 0);

  temple::random::shuffle(indices, engine);

  const unsigned M = boost::num_vertices(_graph);
  std::vector<double> distances (M);
  using ColorMapType = boost::two_bit_color_map<>;
  ColorMapType color_map {M};
  std::vector<VertexDescriptor> predecessors (M);

  std::vector<AtomIndex>::const_iterator separator;

  if(partiality == Partiality::FourAtom) {
    separator = indices.cbegin() + std::min(N, 4U);
  } else if(partiality == Partiality::TenPercent) {
    /* Advance the separator at most by N positions (guards against advancing
     * past-the-end), and at least by four atoms (guards against situations
     * where 4 <= N < 40 and 10% would yield less smoothing than FourAtom)
     *
     * Not equivalent to std::clamp(4u, cast<u>(0.1 * N), N) if N < 4!
     */
    separator = indices.cbegin() + std::min(N, std::max(4U, static_cast<unsigned>(0.1 * N)));
  } else { // All
    separator = indices.cend();
  }

  for(auto iter = indices.cbegin(); iter != separator; ++iter) {
    const AtomIndex a = *iter;

    std::vector<AtomIndex> otherIndices;
    otherIndices.reserve(N - 1);

    // Avoid already-chosen elements
    for(AtomIndex b = 0; b < a; ++b) {
      if(upperTriangle(b, a) == 0) {
        otherIndices.push_back(b);
      }
    }

    for(AtomIndex b = a + 1; b < N; ++b) {
      if(upperTriangle(a, b) == 0) {
        otherIndices.push_back(b);
      }
    }

    temple::random::shuffle(otherIndices, engine);

    // Again through N - 1 indices: N²
    for(const auto& b : otherIndices) {

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
        color_map.data.get() + (color_map.n + ColorMapType::elements_per_char - 1)
          / ColorMapType::elements_per_char,
        0
      );

#ifdef MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
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

      double lower = -distances.at(right(b));
      double upper = distances.at(left(b));

      if(upper < lower) {
        explainContradictionPaths(a, b, predecessors, distances);
        return DgError::GraphImpossible;
      }

      /* Shortest distance from left a vertex to right b vertex is lower bound (negative)
       * Shortest distance from left a vertex to left b vertex is upper bound
       */
      double tightenedBound = temple::random::getSingle<double>(
        lower,
        upper,
        engine
      );

      upperTriangle(
        std::min(a, b),
        std::max(a, b)
      ) = tightenedBound;

      // Modify the graph accordingly
      _updateGraphWithFixedDistance(a, b, tightenedBound);
    }
  }

  for(auto iter = separator; iter != indices.cend(); ++iter) {
    const AtomIndex a = *iter;

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
      color_map.data.get() + (color_map.n + ColorMapType::elements_per_char - 1)
        / ColorMapType::elements_per_char,
      0
    );

#ifdef MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
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

    for(AtomIndex b = 0; b < N; ++b) {
      if(a == b || upperTriangle(std::min(a, b), std::max(a, b)) > 0) {
        continue;
      }

      double presumedLower = -distances.at(right(b));
      double presumedUpper = distances.at(left(b));

      /* Shortest distance from left a vertex to right b vertex is lower bound (negative)
       * Shortest distance from left a vertex to left b vertex is upper bound
       */
      double tightenedBound = temple::random::getSingle<double>(
        std::min(presumedLower, presumedUpper),
        std::max(presumedLower, presumedUpper),
        engine
      );

      upperTriangle(
        std::min(a, b),
        std::max(a, b)
      ) = tightenedBound;

      // Modify the graph accordingly
      _updateGraphWithFixedDistance(a, b, tightenedBound);
    }
  }

  return distancesMatrix;
}

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine
