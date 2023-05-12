/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/DistanceGeometry/ExplicitBoundsGraph.h"

#include "boost/graph/two_bit_color_map.hpp"

#include "Molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "Molassembler/DistanceGeometry/DistanceGeometry.h"
#include "Molassembler/DistanceGeometry/Error.h"
#include "Molassembler/Log.h"
#include "Molassembler/Modeling/AtomInfo.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Options.h"
#include "Molassembler/Graph/PrivateGraph.h"

#include "Molassembler/Temple/Random.h"

#ifdef MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
#include "Molassembler/DistanceGeometry/Gor1.h"
#else
#include "Molassembler/Graph/Gor1.h"
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
namespace Molassembler {
namespace DistanceGeometry {

ExplicitBoundsGraph::ExplicitBoundsGraph(
  const PrivateGraph& inner,
  const BoundsMatrix& bounds
) : graph_ {2 * inner.V()},
    inner_ {inner}
{
  const AtomIndex N = inner.V();

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

        boost::add_edge(left(a), right(b), -vdwLowerBound, graph_);
        boost::add_edge(left(b), right(a), -vdwLowerBound, graph_);
      } else {
        /* Add explicit edges for this a-b pair */
        // Bidirectional edge in left graph with upper weight
        boost::add_edge(left(a), left(b), upperBound, graph_);
        boost::add_edge(left(b), left(a), upperBound, graph_);

        // Bidirectional edge in right graph with upper weight
        boost::add_edge(right(a), right(b), upperBound, graph_);
        boost::add_edge(right(b), right(a), upperBound, graph_);

        // Forward edge from left to right graph with negative lower bound weight
        boost::add_edge(left(a), right(b), -lowerBound, graph_);
        boost::add_edge(left(b), right(a), -lowerBound, graph_);
      }
    }
  }

  // Determine the two heaviest element types in the molecule, O(N)
  heaviestAtoms_ = {{Utils::ElementType::H, Utils::ElementType::H}};
  for(AtomIndex i = 0; i < N; ++i) {
    auto elementType = inner.elementType(i);
    if(
      Utils::ElementInfo::Z(elementType)
      > Utils::ElementInfo::Z(heaviestAtoms_.back())
    ) {
      heaviestAtoms_.back() = elementType;

      if(
        Utils::ElementInfo::Z(heaviestAtoms_.back())
        > Utils::ElementInfo::Z(heaviestAtoms_.front())
      ) {
        std::swap(heaviestAtoms_.front(), heaviestAtoms_.back());
      }
    }
  }
}

ExplicitBoundsGraph::ExplicitBoundsGraph(
  const PrivateGraph& inner,
  const DistanceBoundsMatrix& bounds
) : graph_ {2 * inner.V()},
    inner_ {inner}
{
  const VertexDescriptor N = inner.V();
  for(VertexDescriptor a = 0; a < N; ++a) {
    for(VertexDescriptor b = a + 1; b < N; ++b) {
      const double lower = bounds.lowerBound(a, b);
      const double upper = bounds.upperBound(a, b);

      if(lower != DistanceBoundsMatrix::defaultLower) {
        // Forward edge from left to right graph with negative lower bound weight
        boost::add_edge(left(a), right(b), -lower, graph_);
        boost::add_edge(left(b), right(a), -lower, graph_);
      } else {
        const double vdwLowerBound = (
          AtomInfo::vdwRadius(inner_.elementType(a))
          + AtomInfo::vdwRadius(inner_.elementType(b))
        );

        // Implicit lower bound on distance between the vertices
        boost::add_edge(left(a), right(b), -vdwLowerBound, graph_);
        boost::add_edge(left(b), right(a), -vdwLowerBound, graph_);
      }

      if(upper != DistanceBoundsMatrix::defaultUpper) {
        // Bidirectional edge in left graph with upper weight
        boost::add_edge(left(a), left(b), upper, graph_);
        boost::add_edge(left(b), left(a), upper, graph_);

        // Bidirectional edge in right graph with upper weight
        boost::add_edge(right(a), right(b), upper, graph_);
        boost::add_edge(right(b), right(a), upper, graph_);
      }
    }
  }

  // Determine the two heaviest element types in the molecule, O(N)
  heaviestAtoms_ = {{Utils::ElementType::H, Utils::ElementType::H}};
  for(AtomIndex i = 0; i < N; ++i) {
    auto elementType = inner.elementType(i);
    if(
      Utils::ElementInfo::Z(elementType)
      > Utils::ElementInfo::Z(heaviestAtoms_.back())
    ) {
      heaviestAtoms_.back() = elementType;

      if(
        Utils::ElementInfo::Z(heaviestAtoms_.back())
        > Utils::ElementInfo::Z(heaviestAtoms_.front())
      ) {
        std::swap(heaviestAtoms_.front(), heaviestAtoms_.back());
      }
    }
  }
}

void ExplicitBoundsGraph::addBound(
  const VertexDescriptor a,
  const VertexDescriptor b,
  const ValueBounds& bound
) {
  // Bidirectional edge in left graph with upper weight
  boost::add_edge(left(a), left(b), bound.upper, graph_);
  boost::add_edge(left(b), left(a), bound.upper, graph_);

  // Bidirectional edge in right graph with upper weight
  boost::add_edge(right(a), right(b), bound.upper, graph_);
  boost::add_edge(right(b), right(a), bound.upper, graph_);

  // Forward edge from left to right graph with negative lower bound weight
  boost::add_edge(left(a), right(b), -bound.lower, graph_);
  boost::add_edge(left(b), right(a), -bound.lower, graph_);
}

void ExplicitBoundsGraph::explainContradictionPaths(
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

void ExplicitBoundsGraph::updateOrAddEdge_(
  const VertexDescriptor i,
  const VertexDescriptor j,
  const double edgeWeight
) {
  auto edgeSearchPair = boost::edge(i, j, graph_);
  if(edgeSearchPair.second) {
    boost::get(boost::edge_weight, graph_, edgeSearchPair.first) = edgeWeight;
  } else {
    boost::add_edge(i, j, edgeWeight, graph_);
  }
}

void ExplicitBoundsGraph::updateGraphWithFixedDistance_(
  const VertexDescriptor a,
  const VertexDescriptor b,
  const double fixedDistance
) {
  updateOrAddEdge_(left(a), left(b), fixedDistance);
  updateOrAddEdge_(left(b), left(a), fixedDistance);

  updateOrAddEdge_(right(a), right(b), fixedDistance);
  updateOrAddEdge_(right(b), right(a), fixedDistance);

  updateOrAddEdge_(left(a), right(b), -fixedDistance);
  updateOrAddEdge_(left(b), right(a), -fixedDistance);
}

double ExplicitBoundsGraph::lowerBound(
  const VertexDescriptor a,
  const VertexDescriptor b
) const {
  auto edgeSearchPair = boost::edge(left(a), right(b), graph_);

  assert(edgeSearchPair.second);

  // The graph contains the lower bound negated
  return -boost::get(boost::edge_weight, graph_, edgeSearchPair.first);
}

double ExplicitBoundsGraph::upperBound(
  const VertexDescriptor a,
  const VertexDescriptor b
) const {
  auto edgeSearchPair = boost::edge(left(a), left(b), graph_);

  assert(edgeSearchPair.second);

  return boost::get(boost::edge_weight, graph_, edgeSearchPair.first);
}

double ExplicitBoundsGraph::maximalImplicitLowerBound(const VertexDescriptor i) const {
  assert(isLeft(i));
  AtomIndex a = i / 2;
  Utils::ElementType elementType = inner_.elementType(a);

  if(elementType == heaviestAtoms_.front()) {
    return AtomInfo::vdwRadius(
      heaviestAtoms_.back()
    ) + AtomInfo::vdwRadius(elementType);
  }

  return AtomInfo::vdwRadius(
    heaviestAtoms_.front()
  ) + AtomInfo::vdwRadius(elementType);
}

const ExplicitBoundsGraph::GraphType& ExplicitBoundsGraph::graph() const {
  return graph_;
}

Result<Eigen::MatrixXd> ExplicitBoundsGraph::makeDistanceBounds() const noexcept {
  unsigned N = inner_.V();

  Eigen::MatrixXd bounds;
  bounds.resize(N, N);
  bounds.setZero();

  unsigned M = boost::num_vertices(graph_);
  std::vector<double> distances (M);
  std::vector<VertexDescriptor> predecessors (M);
  using ColorMapType = boost::two_bit_color_map<>;
  ColorMapType color_map {M};

  for(AtomIndex a = 0; a < N - 1; ++a) {
    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      boost::get(boost::vertex_index, graph_)
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      boost::get(boost::vertex_index, graph_)
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
      graph_,
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

Result<Eigen::MatrixXd> ExplicitBoundsGraph::makeDistanceMatrix(Random::Engine& engine) noexcept {
  return makeDistanceMatrix(engine, Partiality::All);
}

Result<Eigen::MatrixXd> ExplicitBoundsGraph::makeDistanceMatrix(Random::Engine& engine, Partiality partiality) noexcept {
  const unsigned N = inner_.V();

  Eigen::MatrixXd distancesMatrix;
  distancesMatrix.resize(N, N);
  distancesMatrix.setZero();

  auto upperTriangle = distancesMatrix.triangularView<Eigen::StrictlyUpper>();

  std::vector<AtomIndex> indices (N);

  std::iota(std::begin(indices), std::end(indices), 0);

  Temple::Random::shuffle(indices, engine);

  const unsigned M = boost::num_vertices(graph_);
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

    Temple::Random::shuffle(otherIndices, engine);

    // Again through N - 1 indices: N²
    for(const auto& b : otherIndices) {

      auto predecessor_map = boost::make_iterator_property_map(
        predecessors.begin(),
        boost::get(boost::vertex_index, graph_)
      );

      auto distance_map = boost::make_iterator_property_map(
        distances.begin(),
        boost::get(boost::vertex_index, graph_)
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
        graph_,
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
      double tightenedBound = Temple::Random::getSingle<double>(
        lower,
        upper,
        engine
      );

      upperTriangle(
        std::min(a, b),
        std::max(a, b)
      ) = tightenedBound;

      // Modify the graph accordingly
      updateGraphWithFixedDistance_(a, b, tightenedBound);
    }
  }

  for(auto iter = separator; iter != indices.cend(); ++iter) {
    const AtomIndex a = *iter;

    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      boost::get(boost::vertex_index, graph_)
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      boost::get(boost::vertex_index, graph_)
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
      graph_,
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
      double tightenedBound = Temple::Random::getSingle<double>(
        std::min(presumedLower, presumedUpper),
        std::max(presumedLower, presumedUpper),
        engine
      );

      upperTriangle(
        std::min(a, b),
        std::max(a, b)
      ) = tightenedBound;

      // Modify the graph accordingly
      updateGraphWithFixedDistance_(a, b, tightenedBound);
    }
  }

  return distancesMatrix;
}

} // namespace DistanceGeometry
} // namespace Molassembler
} // namespace Scine
