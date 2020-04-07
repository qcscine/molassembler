/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "molassembler/DistanceGeometry/ImplicitBoundsGraphBoost.h"
#include "boost/graph/two_bit_color_map.hpp"

#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/DistanceGeometry/Error.h"
#include "molassembler/Graph/PrivateGraph.h"
#include "molassembler/Log.h"
#include "molassembler/Modeling/AtomInfo.h"
#include "molassembler/Options.h"

#include "molassembler/Temple/Random.h"

#include "Utils/Geometry/ElementInfo.h"

#ifdef MOLASSEMBLER_IMPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
#include "molassembler/DistanceGeometry/Gor1.h"
#else
#include "molassembler/Graph/Gor1.h"
#endif


namespace Scine {
namespace molassembler {
namespace distance_geometry {

/* Class Implementation */

void ImplicitBoundsGraph::explainContradictionPaths_(
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
    logRef << "Encountered contradiction in triangle ineqaulity limits calculation.\n";
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

ImplicitBoundsGraph::ImplicitBoundsGraph(
  const PrivateGraph& inner,
  BoundsMatrix bounds
) : innerGraphPtr_(&inner), distances_(std::move(bounds)) {
  // Determine the two heaviest element types in the molecule, O(N)
  heaviestAtoms_ = {{Utils::ElementType::H, Utils::ElementType::H}};
  const VertexDescriptor N = inner.N();
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

ImplicitBoundsGraph::VertexDescriptor ImplicitBoundsGraph::num_vertices() const {
  return 2 * distances_.outerSize();
}

ImplicitBoundsGraph::VertexDescriptor ImplicitBoundsGraph::num_edges() const {
  /* Every entry in the upper triangle indicates a bound
   * For every upper bound, there are 4 edges
   * For every lower bound, there are 2 edges
   * If there is no entry in the upper triangle, there is still an implicit
   * lower bound.
   */
  unsigned N = distances_.outerSize();
  unsigned count = 0;

  for(unsigned i = 0; i < N - 1; ++i) {
    for(unsigned j = i + 1; j < N; ++j) {
      if(distances_(i, j) > 0) {
        ++count;
      }
    }
  }

  return 4 * count + N * (N - 1);
}

std::pair<ImplicitBoundsGraph::EdgeDescriptor, bool> ImplicitBoundsGraph::edge(const VertexDescriptor i, const VertexDescriptor j) const {
  const auto a = internal(i);
  const auto b = internal(j);

  // Right-to-left edges never exist
  if(!isLeft(i) && isLeft(j)) {
    return {
      EdgeDescriptor {},
      false
    };
  }

  // Left-to-right edges always exist, whether explicit or implicit does not matter
  if(isLeft(i) && !isLeft(j) && a != b) {
    return {
      EdgeDescriptor {i, j},
      true
    };
  }

  unsigned N = distances_.outerSize();

  if(a < N && b < N && distances_(a, b) != 0) {
    return {
      EdgeDescriptor {i, j},
      true
    };
  }

  return {
    EdgeDescriptor {},
    false
  };
}

bool ImplicitBoundsGraph::hasExplicit(const EdgeDescriptor& edge) const {
  assert(isLeft(edge.first) && !isLeft(edge.second));

  return distances_(internal(edge.first), internal(edge.second)) != 0;
}

outcome::result<Eigen::MatrixXd> ImplicitBoundsGraph::makeDistanceBounds() const noexcept {
  Eigen::MatrixXd bounds;

  unsigned N = distances_.outerSize();
  bounds.resize(N, N);
  bounds.setZero();

  unsigned M = num_vertices();
  std::vector<double> distances (M);
  std::vector<VertexDescriptor> predecessors (M);
  using ColorMapType = boost::two_bit_color_map<>;
  ColorMapType color_map {M};

  for(VertexDescriptor a = 0; a < N; ++a) {
    // Perform a single shortest paths calculation for a

    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      VertexIndexMap()
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      VertexIndexMap()
    );

    // re-fill color map with white
    std::fill(
      color_map.data.get(),
      color_map.data.get() + (color_map.n + ColorMapType::elements_per_char - 1)
        / ColorMapType::elements_per_char,
      0
    );

#ifdef MOLASSEMBLER_IMPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
    boost::gor1_ig_shortest_paths(
      *this,
      VertexDescriptor {left(a)},
      predecessor_map,
      color_map,
      distance_map
    );
#else
    boost::gor1_simplified_shortest_paths(
      *this,
      VertexDescriptor {left(a)},
      predecessor_map,
      color_map,
      distance_map
    );
#endif

    for(VertexDescriptor b = a + 1; b < N; ++b) {
      // a is always smaller than b, hence (a, b) is the upper bound
      bounds(a, b) = distances.at(left(b));
      bounds(b, a) = -distances.at(right(b));

      if(bounds(a, b) < bounds(b, a)) {
        explainContradictionPaths_(a, b, predecessors, distances);
        return DgError::GraphImpossible;
      }
    }
  }

  return bounds;
}

outcome::result<Eigen::MatrixXd> ImplicitBoundsGraph::makeDistanceMatrix(random::Engine& engine) noexcept {
  return makeDistanceMatrix(engine, Partiality::All);
}

outcome::result<Eigen::MatrixXd> ImplicitBoundsGraph::makeDistanceMatrix(random::Engine& engine, Partiality partiality) noexcept {
  const unsigned N = innerGraphPtr_->N();

  std::vector<AtomIndex> indices(N);
  std::iota(
    indices.begin(),
    indices.end(),
    0
  );

  temple::random::shuffle(indices, engine);

  const unsigned M = num_vertices();
  std::vector<double> distances (M);
  std::vector<VertexDescriptor> predecessors (M);

  using ColorMapType = boost::two_bit_color_map<>;
  // Determine triangle inequality limits for pair
  ColorMapType color_map {M};

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

    for(AtomIndex b = 0; b < a; ++b) {
      if(!(distances_(a, b) == distances_(b, a) && distances_(a, b) != 0)) {
        otherIndices.push_back(b);
      }
    }
    for(AtomIndex b = a + 1; b < N; ++b) {
      if(!(distances_(a, b) == distances_(b, a) && distances_(a, b) != 0)) {
        otherIndices.push_back(b);
      }
    }

    temple::random::shuffle(otherIndices, engine);

    for(const AtomIndex b : otherIndices) {
      auto predecessor_map = boost::make_iterator_property_map(
        predecessors.begin(),
        VertexIndexMap()
      );

      auto distance_map = boost::make_iterator_property_map(
        distances.begin(),
        VertexIndexMap()
      );

      // re-fill color map with white
      std::fill(
        color_map.data.get(),
        color_map.data.get() + (color_map.n + ColorMapType::elements_per_char - 1)
          / ColorMapType::elements_per_char,
        0
      );

#ifdef MOLASSEMBLER_IMPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
      boost::gor1_ig_shortest_paths(
        *this,
        VertexDescriptor {left(a)},
        predecessor_map,
        color_map,
        distance_map
      );

      double presumedLower = -distances.at(right(b));
      double presumedUpper = distances.at(left(b));

      double fixedDistance = temple::random::getSingle<double>(
        std::min(presumedLower, presumedUpper),
        std::max(presumedLower, presumedUpper),
        engine
      );
#else
      boost::gor1_simplified_shortest_paths(
        *this,
        VertexDescriptor {left(a)},
        predecessor_map,
        color_map,
        distance_map
      );

      if(distances.at(left(b)) < -distances.at(right(b))) {
        explainContradictionPaths_(a, b, predecessors, distances);
        return DgError::GraphImpossible;
      }

      // Pick fixed distance
      double fixedDistance = temple::random::getSingle<double>(
        -distances.at(right(b)),
        distances.at(left(b)),
        engine
      );

#endif

      // Record in distances matrix
      distances_(a, b) = fixedDistance;
      distances_(b, a) = fixedDistance;
    }
  }

  for(auto iter = separator; iter != indices.cend(); ++iter) {
    const AtomIndex a = *iter;

    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      VertexIndexMap()
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      VertexIndexMap()
    );

    // re-fill color map with white
    std::fill(
      color_map.data.get(),
      color_map.data.get() + (color_map.n + ColorMapType::elements_per_char - 1)
        / ColorMapType::elements_per_char,
      0
    );

#ifdef MOLASSEMBLER_IMPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
    boost::gor1_ig_shortest_paths(
      *this,
      VertexDescriptor {left(a)},
      predecessor_map,
      color_map,
      distance_map
    );
#else
    boost::gor1_simplified_shortest_paths(
      *this,
      VertexDescriptor {left(a)},
      predecessor_map,
      color_map,
      distance_map
    );
#endif

    for(AtomIndex b = 0; b < N; ++b) {
      if(
        a == b
        || (
          distances_(a, b) == distances_(b, a)
          && distances_(a, b) != 0
        )
      ) {
        // Skip on-diagonal and already-chosen entries
        continue;
      }

      double presumedLower = -distances.at(right(b));
      double presumedUpper = distances.at(left(b));

      // Pick fixed distance
      double fixedDistance = temple::random::getSingle<double>(
        std::min(presumedLower, presumedUpper),
        std::max(presumedLower, presumedUpper),
        engine
      );

      // Record in distances matrix
      distances_(a, b) = fixedDistance;
      distances_(b, a) = fixedDistance;
    }
  }

  return distances_;
}

ImplicitBoundsGraph::VertexDescriptor ImplicitBoundsGraph::out_degree(VertexDescriptor i) const {
  unsigned count = 0;
  unsigned N = distances_.outerSize();
  VertexDescriptor a = internal(i);
  for(VertexDescriptor b = 0; b < N; ++b) {
    if(distances_(a, b) != 0) {
      ++count;
    }
  }

  if(isLeft(i)) {
    return count + N - 1;
  }

  return count;
}

double& ImplicitBoundsGraph::lowerBound(const VertexDescriptor a, const VertexDescriptor b) {
  return distances_(
    std::max(a, b),
    std::min(a, b)
  );
}

double& ImplicitBoundsGraph::upperBound(const VertexDescriptor a, const VertexDescriptor b) {
  return distances_(
    std::min(a, b),
    std::max(a, b)
  );
}

double ImplicitBoundsGraph::lowerBound(const VertexDescriptor a, const VertexDescriptor b) const {
  return distances_(
    std::max(a, b),
    std::min(a, b)
  );
}

double ImplicitBoundsGraph::upperBound(const VertexDescriptor a, const VertexDescriptor b) const {
  return distances_(
    std::min(a, b),
    std::max(a, b)
  );
}

double ImplicitBoundsGraph::maximalImplicitLowerBound(const VertexDescriptor i) const {
  assert(isLeft(i));
  auto a = internal(i);
  auto elementType = innerGraphPtr_->elementType(a);

  if(elementType == heaviestAtoms_.front()) {
    return atom_info::vdwRadius(
      heaviestAtoms_.back()
    ) + atom_info::vdwRadius(elementType);
  }

  return atom_info::vdwRadius(
    heaviestAtoms_.front()
  ) + atom_info::vdwRadius(elementType);
}

/* Nested classes */
/* Edge weight property map */
ImplicitBoundsGraph::EdgeWeightMap::EdgeWeightMap(const ImplicitBoundsGraph& base)
  : basePtr_(&base)
{}

double ImplicitBoundsGraph::EdgeWeightMap::operator [] (const EdgeDescriptor& e) const {
  auto a = internal(e.first);
  auto b = internal(e.second);

  if(isLeft(e.first) && !isLeft(e.second)) {
    double explicitValue = basePtr_->lowerBound(a, b);

    if(explicitValue != 0.0) {
      return -explicitValue;
    }

    return -(
      atom_info::vdwRadius(
        basePtr_->innerGraphPtr_->elementType(a)
      ) + atom_info::vdwRadius(
        basePtr_->innerGraphPtr_->elementType(b)
      )
    );
  }

  return basePtr_->upperBound(a, b);
}

double ImplicitBoundsGraph::EdgeWeightMap::operator () (const EdgeDescriptor& e) const {
  return this->operator[](e);
}

ImplicitBoundsGraph::EdgeWeightMap ImplicitBoundsGraph::getEdgeWeightPropertyMap() const {
  return {*this};
}

/* Vertex iterator */
ImplicitBoundsGraph::vertex_iterator::vertex_iterator() = default;
ImplicitBoundsGraph::vertex_iterator::vertex_iterator(ImplicitBoundsGraph::VertexDescriptor i) : index(i) {}
ImplicitBoundsGraph::vertex_iterator::vertex_iterator(const ImplicitBoundsGraph::vertex_iterator& other) = default;
ImplicitBoundsGraph::vertex_iterator::vertex_iterator(ImplicitBoundsGraph::vertex_iterator&& other) noexcept = default;
ImplicitBoundsGraph::vertex_iterator& ImplicitBoundsGraph::vertex_iterator::operator = (const ImplicitBoundsGraph::vertex_iterator& other) = default;
ImplicitBoundsGraph::vertex_iterator& ImplicitBoundsGraph::vertex_iterator::operator = (ImplicitBoundsGraph::vertex_iterator&& other) noexcept = default;

bool ImplicitBoundsGraph::vertex_iterator::operator == (const ImplicitBoundsGraph::vertex_iterator& other) const {
  return index == other.index;
}

bool ImplicitBoundsGraph::vertex_iterator::operator != (const ImplicitBoundsGraph::vertex_iterator& other) const {
  return index != other.index;
}

ImplicitBoundsGraph::vertex_iterator& ImplicitBoundsGraph::vertex_iterator::operator ++ () {
  ++index;
  return *this;
}

ImplicitBoundsGraph::vertex_iterator ImplicitBoundsGraph::vertex_iterator::operator ++ (int) {
  vertex_iterator copy = *this;
  ++index;
  return copy;
}

ImplicitBoundsGraph::vertex_iterator& ImplicitBoundsGraph::vertex_iterator::operator -- () {
  --index;
  return *this;
}

ImplicitBoundsGraph::vertex_iterator ImplicitBoundsGraph::vertex_iterator::operator -- (int) {
  vertex_iterator copy = *this;
  --index;
  return copy;
}

ImplicitBoundsGraph::vertex_iterator& ImplicitBoundsGraph::vertex_iterator::operator + (unsigned i) {
  index += i;
  return *this;
}

ImplicitBoundsGraph::vertex_iterator& ImplicitBoundsGraph::vertex_iterator::operator - (unsigned i) {
  index -= i;
  return *this;
}

ImplicitBoundsGraph::VertexDescriptor ImplicitBoundsGraph::vertex_iterator::operator * () const {
  return index;
}

int ImplicitBoundsGraph::vertex_iterator::operator - (const vertex_iterator& other) const {
  return static_cast<int>(index) - static_cast<int>(other.index);
}

ImplicitBoundsGraph::vertex_iterator ImplicitBoundsGraph::vbegin() const {
  return {};
}

ImplicitBoundsGraph::vertex_iterator ImplicitBoundsGraph::vend() const {
  return {num_vertices()};
}

/* Edge iterator */
ImplicitBoundsGraph::edge_iterator::edge_iterator() = default;
ImplicitBoundsGraph::edge_iterator::edge_iterator(
  const ImplicitBoundsGraph& base,
  VertexDescriptor i
) : basePtr_ {&base},
    i_ {i}
{
  auto a = internal(i_);

  b_ = 0;
  if(a == 0) {
    ++b_;
  }

  crossGroup_ = false;
  unsigned N = basePtr_->distances_.outerSize();

  if(a < N) {
    while(b_ < N && basePtr_->distances_(a, b_) == 0.0) {
      ++b_;
    }
  }
}

ImplicitBoundsGraph::edge_iterator::edge_iterator(const ImplicitBoundsGraph::edge_iterator& other) = default;
ImplicitBoundsGraph::edge_iterator::edge_iterator(ImplicitBoundsGraph::edge_iterator&& other) noexcept = default;
ImplicitBoundsGraph::edge_iterator& ImplicitBoundsGraph::edge_iterator::operator = (const ImplicitBoundsGraph::edge_iterator& other) = default;
ImplicitBoundsGraph::edge_iterator& ImplicitBoundsGraph::edge_iterator::operator = (ImplicitBoundsGraph::edge_iterator&& other) noexcept = default;

ImplicitBoundsGraph::edge_iterator& ImplicitBoundsGraph::edge_iterator::operator ++ () {
  increment_();
  return *this;
}

ImplicitBoundsGraph::edge_iterator ImplicitBoundsGraph::edge_iterator::operator ++ (int) {
  edge_iterator copy = *this;
  ++(*this);
  return copy;
}

bool ImplicitBoundsGraph::edge_iterator::operator == (const edge_iterator& other) const {
  return (
    crossGroup_ == other.crossGroup_
    && b_ == other.b_
    && i_ == other.i_
    && basePtr_ == other.basePtr_
  );
}

bool ImplicitBoundsGraph::edge_iterator::operator != (const edge_iterator& other) const {
  return !(*this == other);
}

double ImplicitBoundsGraph::edge_iterator::weight() const {
  auto a = internal(i_);
  if(crossGroup_) {
    double data = basePtr_->lowerBound(a, b_);

    if(data != 0.0) {
      return -data;
    }

    return -(
      atom_info::vdwRadius(
        basePtr_->innerGraphPtr_->elementType(a)
      ) + atom_info::vdwRadius(
        basePtr_->innerGraphPtr_->elementType(b_)
      )
    );
  }

  return basePtr_->upperBound(a, b_);
}

ImplicitBoundsGraph::VertexDescriptor ImplicitBoundsGraph::edge_iterator::target() const {
  if(crossGroup_) {
    return right(b_);
  }

  if(isLeft(i_)) {
    return left(b_);
  }

  return right(b_);
}

ImplicitBoundsGraph::EdgeDescriptor ImplicitBoundsGraph::edge_iterator::operator * () const {
  return {
    i_,
    target()
  };
}

void ImplicitBoundsGraph::edge_iterator::increment_() {
  unsigned N = basePtr_->distances_.outerSize();
  auto a = internal(i_);

  if(!crossGroup_) {
    ++b_;
    // Find the next explicit information
    while(b_ < N && basePtr_->distances_(a, b_) == 0.0) {
      ++b_;
    }

    if(b_ == N) {
      if(isLeft(i_)) {
        // Roll over to implicits only if is_ is left
        crossGroup_ = true;
        b_ = 0;
        if(b_ == a) {
          ++b_;
        }
      } else {
        // Right vertices have no implicits
        ++i_;
        b_ = 0;

        a = internal(i_);

        if(b_ == a) {
          ++b_;
        }

        if(a < N) {
          // Search for next explicit
          while(b_ < N && basePtr_->distances_(a, b_) == 0.0) {
            ++b_;
          }
        }
      }
    }
  } else {
    ++b_;

    if(b_ == a) {
      ++b_;
    }

    if(b_ == N) {
      // Rollover to next explicits
      crossGroup_ = false;
      b_ = 0;
      ++i_;

      a = internal(i_);

      if(b_ == a) {
        ++b_;
      }

      while(b_ < N && basePtr_->distances_(a, b_) == 0.0) {
        ++b_;
      }
    }
  }
}

std::string ImplicitBoundsGraph::edge_iterator::state() const {
  using namespace std::string_literals;

  return std::to_string(i_) + " "s
    + std::to_string(b_) + " "s
    + std::to_string(static_cast<int>(crossGroup_));
}

ImplicitBoundsGraph::edge_iterator ImplicitBoundsGraph::ebegin() const {
  return {
    *this,
    0
  };
}

ImplicitBoundsGraph::edge_iterator ImplicitBoundsGraph::eend() const {
  return {
    *this,
    num_vertices()
  };
}


/* Out edge iterator */
ImplicitBoundsGraph::edge_iterator ImplicitBoundsGraph::obegin(VertexDescriptor i) const {
  return {
    *this,
    i
  };
}

ImplicitBoundsGraph::edge_iterator ImplicitBoundsGraph::oend(VertexDescriptor i) const {
  return {
    *this,
    i + 1
  };
}

ImplicitBoundsGraph::in_group_edge_iterator::in_group_edge_iterator() = default;
ImplicitBoundsGraph::in_group_edge_iterator::in_group_edge_iterator(
  const ImplicitBoundsGraph::in_group_edge_iterator& other
) = default;
ImplicitBoundsGraph::in_group_edge_iterator::in_group_edge_iterator(
  ImplicitBoundsGraph::in_group_edge_iterator&& other
) noexcept = default;
ImplicitBoundsGraph::in_group_edge_iterator& ImplicitBoundsGraph::in_group_edge_iterator::operator = (
  ImplicitBoundsGraph::in_group_edge_iterator&& other
) noexcept = default;
ImplicitBoundsGraph::in_group_edge_iterator& ImplicitBoundsGraph::in_group_edge_iterator::operator = (
  const ImplicitBoundsGraph::in_group_edge_iterator& other
) = default;

ImplicitBoundsGraph::in_group_edge_iterator::in_group_edge_iterator(
  const ImplicitBoundsGraph& base,
  const VertexDescriptor i
) : basePtr_{&base},
    i_ {i},
    b_ {0},
    isLeft_ {isLeft(i)}
{
  auto a = internal(i_);

  if(a == 0) {
    ++b_;
  }

  unsigned N = basePtr_->distances_.outerSize();
  while(b_ < N && basePtr_->distances_(a, b_) == 0.0) {
    ++b_;
  }
}

ImplicitBoundsGraph::in_group_edge_iterator::in_group_edge_iterator(
  const ImplicitBoundsGraph& base,
  const VertexDescriptor i,
  bool /* tag */
) : basePtr_{&base},
    i_ {i},
    b_ {static_cast<VertexDescriptor>(base.distances_.outerSize())},
    isLeft_ {isLeft(i)}
{
  if(internal(i_) == b_) {
    ++b_;
  }
}

} // namespace distance_geometry
} // namespace molassembler
} // namespace Scine
