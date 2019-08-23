/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Space-efficient shortest paths graph for triangle inequalities
 *
 * Declaration of a graph class that aids in the fast determination of triangle
 * inequality limit distance bounds and the generation of a distance matrix from
 * a list of atom-index pairwise distance bounds.
 */

#ifndef INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_BOUNDS_GRAPH_H
#define INCLUDE_DISTANCE_GEOMETRY_IMPLICIT_BOUNDS_GRAPH_H

// #define MOLASSEMBLER_IMPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM

#include "boost_outcome/outcome.hpp"
#include "boost/property_map/property_map.hpp"
#include "Eigen/Core"
#include "Utils/Geometry/ElementTypes.h"

#include "molassembler/DistanceGeometry/ValueBounds.h"

#include <array>
#include <map>
#include <tuple>

namespace Scine {

namespace molassembler {

namespace random {
class Engine;
} // namespace random

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

// Forward-declare InnerGraph
class InnerGraph;

namespace DistanceGeometry {

// Forward-declarations
enum class Partiality;
class DistanceBoundsMatrix;

/*!
 * @brief Simulates a graph from which triangle inequality bounds can be calculated by shortest-paths
 *
 * Based off of pairwise bounds collected from SpatialModel, this class
 * simulates a graph in which one-to-all triangle inequality bounds can be
 * calculated by a single-source shortest paths calculation.
 *
 * The general construction of the graph is loosely described in:
 * Havel, Timothy F. "Distance geometry: Theory, algorithms, and chemical
 * applications." Encyclopedia of Computational Chemistry 120 (1998): 723-742.
 *
 * In essence, for every atomic index a ∈ [0, N), the graph contains two
 * vertices, left(a) and right(a). For every pair of atomic indices a, b that
 * have a distance bound, i.e. a lower and upper bound on their spatial distance
 * within the current metric space, there are six directed edges in the graph:
 * - Two edges (one bidirectional edge) between left(a) and left(b) with the
 *   edge weight of the upper bound
 * - Two edges between right(a) and right(b) with the edge weight of
 *   the upper bound
 * - A unidirectional edge from left(a) to right(b) with the edge weight of
 *   the lower bound, negated
 * - A unidirectional edge from left(b) to right(a) with the edge weight of
 *   the lower bound, negated
 *
 * If there is no explicit bounding information present for the pair a, b, then
 * there are still two edges in the graph:
 * - A unidirectional edge from left(a) to right(b) with the edge weight of
 *   the respective atom types' van-der-Walls radii summed
 * - A unidirectional edge from left(b) to right(a) with the edge weight of
 *   the respective atom types' van-der-Walls radii summed
 *
 * There are never edges from the right part of the graph to the left.
 *
 * If a distance bound is narrowed to a fixed distance, this is represented as
 * raising the lower bound to the fixed distance and lowering the upper bound to
 * the fixed distance.
 *
 * A lot of the implemented functions are not strictly necessary for the
 * specialized shortest-paths algorithm that is used in makeDistanceMatrix(),
 * which is the main purpose of the class, but exist to allow this graph to
 * integrate with boost::graph in order to easily benchmark it with other
 * algorithms like Bellman-Ford and compare it with ExplicitGraph, which uses a
 * fully-explicit boost::graph class as its underlying datastructure.
 *
 * Why is the implicit nature of the bonds on the basis of an underlying matrix
 * so much faster than ExplicitGraph when it comes to generating distance
 * matrices?
 * - The GOR and Bellman-Ford shortest-paths algorithms have a complexity of
 *   O(V * E), where V is the number of vertices and E the number of edges.
 *   Since for every left vertex V, there are V-1 edges to the right side
 *   (whether explicit or implicit does not matter), E is on the order of V²,
 *   so the entire shortest-paths algorithm is O(V³). Under certain conditions
 *   when on a left vertex during the shortest-paths calculation, one can
 *   choose to ignore all edges to the right side. While in an explicit BGL
 *   graph, this would entail getting the target vertex, comparing and
 *   incrementing the iterator, in this class, which merely simulates the
 *   existence of edges, the entire logic of an alternate iterator class that
 *   simulates only in-group (i.e. left-left and right-right) edges can be
 *   significantly reduced. This condition leads to a big reduction in the
 *   number of possible paths and thus heavily accelerates the shortest paths
 *   calculation.
 * - Less memory use is perhaps also a factor. In the generation of distance
 *   bounds, ExplicitGraph requires memory for its underlying BGL graph and edge
 *   weight map, which are heavily redundant, plus an N² double matrix for
 *   storing the resulting distance matrix. ImplicitGraph requires only
 *   the N² double matrix.
 *
 * Notes for reading the implementation:
 * - a, b denote internal indices (molecule [0-N) interval)
 * - i, j denote external indices (since graph simulates 2N vertices)
 * - The functions left, right, isLeft and internal help to interconvert indices
 * - If you encounter an expression like i % 2 == j % 2, this tests whether
 *   both indices are in the same part (left or right) of the graph.
 *
 *
 * Some advice on optimization opportunities:
 *
 * Do not try to apply an algorithm that uses the existing shortest path
 * information from a previous run to calculate the new shortest paths. There
 * are already few shortest paths algorithms well suited to graphs with negative
 * real edge weights. Dynamic algorithms are differentiated into whether they
 * can deal with edge weight changes, vertex additions, vertex removals, etc.
 * A fully dynamic algorithm can deal with any of those (albeit only a singular
 * one at a time). A fully batch dynamic algorithm can deal with an arbitrary
 * number of any type of graph changes. You need one of those, since fixing
 * a distance in ImplicitGraph either entails changing six edge weights
 * or changing one and adding five new edges.
 *
 * These algorithms already have a really hard time getting to better bounds
 * than recalculating the solution from scratch, which is where most effort has
 * been poured into.
 *
 * In the same vein, do not try to decide which parts of the shortest paths
 * tree are volatile after change in order to avoid recalculating the entire
 * solution. The idea here is to avoid recalculating the tree if the next a-b
 * pair is in the stable part of the shortest-paths tree after the prior
 * distance fixture. Creating an oracle to decide this is probably about as hard
 * as creating a batch dynamic SSSP algorithm that can deal with negative edge
 * weights, so: really hard.
 */
class ImplicitGraph {
public:
  using VertexDescriptor = unsigned long;
  using EdgeDescriptor = std::pair<VertexDescriptor, VertexDescriptor>;

private:
  //! Pointer to molecule from which the bounds come from
  const InnerGraph* _innerGraphPtr;

  //! Stores the two heaviest element types
  std::array<Scine::Utils::ElementType, 2> _heaviestAtoms;

  //! Dense adjacency matrix for O(1) access to fixed distances
  Eigen::MatrixXd _distances;

  /* To outer indexing */
  [[gnu::const]] inline static VertexDescriptor left(const VertexDescriptor a) noexcept {
    return 2 * a;
  }

  [[gnu::const]] inline static VertexDescriptor right(const VertexDescriptor a) noexcept {
    return 2 * a + 1;
  }

  [[gnu::const]] inline static VertexDescriptor internal(const VertexDescriptor i) noexcept {
    // Integer division rounds down, which is perfect
    return i / 2;
  }

  static void _explainContradictionPaths(
    VertexDescriptor a,
    VertexDescriptor b,
    const std::vector<VertexDescriptor>& predecessors,
    const std::vector<double>& distances
  );

public:
  using BoundsMatrix = Eigen::MatrixXd;

  ImplicitGraph(
    const InnerGraph& inner,
    BoundsMatrix bounds
  );

  /* Information */
  /*! @brief Returns the number of vertices simulated by the graph, which is 2N
   *
   * Returns the number of vertices simulated by the graph, which is 2 * N
   *
   * @complexity{@math{\Theta(1)}}
   */
  VertexDescriptor num_vertices() const;

  /*! @brief Returns the number of edges currently simulated by the graph
   *
   * Returns the number of edges currently simulated by the graph.
   *
   * @complexity{@math{\Theta(N^2)}}
   */
  VertexDescriptor num_edges() const;

  /*! @brief Fetches an edge descriptor for a speculative edge from i to j
   *
   * Fetches an edge descriptor for the edge from i to j, and whether there is
   * such an edge. If there is no such edge (i.e. the pair's second boolean is
   * false), then using the EdgeDescriptor in other functions results in
   * undefined behavior.
   *
   * @complexity{@math{\Theta(1)}}
   */
  std::pair<EdgeDescriptor, bool> edge(VertexDescriptor i, VertexDescriptor j) const;

  /*! @brief Checks if there is explicit information present for a left-to-right edge
   *
   * @complexity{@math{\Theta(1)}}
   */
  bool hasExplicit(const EdgeDescriptor& edge) const;

  /* To inner indexing */
  /*! @brief Check if a vertex descriptor is part of the left subgraph
   *
   * @complexity{@math{\Theta(1)}}
   */
  [[gnu::const]] inline static bool isLeft(const VertexDescriptor i) noexcept {
    return i % 2 == 0;
  }

  /*! @brief Generates a bare distance bounds matrix by repeated shortest-paths calculations
   *
   * This function generates a matrix representing distance bounds, NOT a
   * DistanceBoundsMatrix.
   *
   * Complexity: O(N * O(shortest paths algorithm))
   * @complexity{@math{\Theta(V^2 \cdot E)}}
   */
  outcome::result<Eigen::MatrixXd> makeDistanceBounds() const noexcept;

  /*! @brief Generates a distance matrix by randomly fixing distances within triangle inequality bounds
   *
   * This function generates a distance matrix. For every pair of atomic indices
   * a, b (traversed at random), a shortest-paths calculation is carried out
   * to determine the triangle inequality limits between a and b. A random
   * distance is chosen between these limits and added to the underlying graph.
   *
   * The generation procedure is destructive, i.e. the same ImplicitGraph
   * cannot be reused to generate multiple distance matrices.
   *
   * @complexity{@math{\Theta(V^3 \cdot E)}. Benchmarked by analysis/BenchmarkGraphAlgorithms}
   *
   * NOTE: This double definition may seem strange, but is necessary to use
   * the forward-declared enum class Partiality correctly.
   */
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(random::Engine& engine) noexcept;

  //!@overload
  outcome::result<Eigen::MatrixXd> makeDistanceMatrix(random::Engine& engine, Partiality partiality) noexcept;

  /*! @brief Returns the source vertex from an edge descriptor
   *
   * @complexity{@math{\Theta(1)}}
   */
  inline VertexDescriptor source(const EdgeDescriptor& e) const {
    return e.first;
  }

  /*! @brief Returns the target vertex from an edge descriptor
   *
   * @complexity{@math{\Theta(1)}}
   */
  inline VertexDescriptor target(const EdgeDescriptor& e) const {
    return e.second;
  }

  /*! @brief Allows const-ref access to the underlying matrix for debug purposes
   *
   * @complexity{@math{\Theta(1)}}
   */
  inline const Eigen::MatrixXd& getMatrix() const {
    return _distances;
  }

  /*! @brief Returns the number of out-edges for a particular vertex
   *
   * @complexity{@math{\Theta(N)}}
   */
  VertexDescriptor out_degree(VertexDescriptor i) const;

  /* These are all O(1) */
  double& lowerBound(VertexDescriptor a, VertexDescriptor b);
  double& upperBound(VertexDescriptor a, VertexDescriptor b);

  double lowerBound(VertexDescriptor a, VertexDescriptor b) const;
  double upperBound(VertexDescriptor a, VertexDescriptor b) const;

  /*! @brief Returns the length of the maximal implicit lower bound outgoing from a left vertex
   *
   * @complexity{@math{\Theta(1)}}
   */
  double maximalImplicitLowerBound(VertexDescriptor i) const;

  //! A helper struct permitting read-access to an edge weight via an edge descriptor
  struct EdgeWeightMap : public boost::put_get_helper<double, EdgeWeightMap> {
    using value_type = double;
    using reference = double;
    using key_type = EdgeDescriptor;

    const ImplicitGraph* _basePtr;

    EdgeWeightMap(const ImplicitGraph& base);
    double operator [] (const EdgeDescriptor& e) const;
    double operator () (const EdgeDescriptor& e) const;
  };

  //! Returns an instance of EdgeWeightMap
  EdgeWeightMap getEdgeWeightPropertyMap() const;

  //! A helper struct to turn vertex descriptors into numeric indices
  struct VertexIndexMap : public boost::put_get_helper<VertexDescriptor, VertexIndexMap> {
    using value_type = VertexDescriptor;
    using key_type = VertexDescriptor;
    using reference = VertexDescriptor;

    inline VertexDescriptor operator [] (const VertexDescriptor v) const { return v; }
    inline VertexDescriptor operator () (const VertexDescriptor v) const { return v; }
  };

  //! A random access iterator through all vertex descriptors of the graph
  class vertex_iterator {
    VertexDescriptor index = 0;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = VertexDescriptor;
    using difference_type = int;
    using pointer = VertexDescriptor*;
    using reference = const VertexDescriptor&;

    vertex_iterator();
    vertex_iterator(VertexDescriptor i);
    vertex_iterator(const vertex_iterator& other);
    vertex_iterator(vertex_iterator&& other) noexcept;
    vertex_iterator& operator = (const vertex_iterator& other);
    vertex_iterator& operator = (vertex_iterator&& other) noexcept;

    bool operator == (const vertex_iterator& other) const;

    bool operator != (const vertex_iterator& other) const;

    vertex_iterator& operator ++ ();

    vertex_iterator operator ++ (int);

    vertex_iterator& operator -- ();

    vertex_iterator operator -- (int);

    vertex_iterator& operator + (unsigned i);

    vertex_iterator& operator - (unsigned i);

    VertexDescriptor operator * () const;

    int operator - (const vertex_iterator& other) const;
  };

  //! Returns a begin-iterator for vertices
  vertex_iterator vbegin () const;

  //! Returns an end-iterator for vertices
  vertex_iterator vend() const;

  //! An iterator to enumerate all edges in a graph.
  class edge_iterator {
    const ImplicitGraph* _basePtr;

    VertexDescriptor _i, _b;
    bool _crossGroup;

    void _increment();

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = EdgeDescriptor;
    using difference_type = int;
    using pointer = EdgeDescriptor*;
    using reference = const EdgeDescriptor&;

    edge_iterator();
    edge_iterator(
      const ImplicitGraph& base,
      VertexDescriptor i
    );
    edge_iterator(const edge_iterator& other);
    edge_iterator(edge_iterator&& other) noexcept;
    edge_iterator& operator = (const edge_iterator& other);
    edge_iterator& operator = (edge_iterator&& other) noexcept;

    edge_iterator& operator ++ ();

    edge_iterator operator ++ (int);

    bool operator == (const edge_iterator& other) const;

    bool operator != (const edge_iterator& other) const;

    std::string state() const;

    EdgeDescriptor operator * () const;

    double weight() const;
    VertexDescriptor target() const;
  };

  //! Returns a begin edge iterator
  edge_iterator ebegin() const;

  //! Returns an end edge iterator
  edge_iterator eend() const;

  //! Returns a begin edge iterator for the out-edges of a specific vertex
  edge_iterator obegin(VertexDescriptor i) const;

  //! Returns an end edge iterator for the out-edges of a specific vertex
  edge_iterator oend(VertexDescriptor i) const;

  //! An iterator to enumerate only edges to the same part of the graph from specific vertices
  class in_group_edge_iterator {
    const ImplicitGraph* _basePtr;
    VertexDescriptor _i;
    VertexDescriptor _b;
    bool _isLeft;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = EdgeDescriptor;
    using difference_type = int;
    using pointer = EdgeDescriptor*;
    using reference = const EdgeDescriptor&;

    // Default constructor
    in_group_edge_iterator();
    // Constructor for the begin iterator
    in_group_edge_iterator(const ImplicitGraph& base, VertexDescriptor i);
    // Constructor for the end iterator
    in_group_edge_iterator(
      const ImplicitGraph& base,
      VertexDescriptor i,
      bool /* tag */
    );
    in_group_edge_iterator(const in_group_edge_iterator& other);
    in_group_edge_iterator(in_group_edge_iterator&& other) noexcept;
    in_group_edge_iterator& operator = (in_group_edge_iterator&& other) noexcept;
    in_group_edge_iterator& operator = (const in_group_edge_iterator& other);

    inline in_group_edge_iterator& operator ++ () {
      /* Find the next position in the matrix with information
       * The fastest way to check if there is information is to just access
       * either the lower or upper bound, and convert the double there to
       * boolean. That boolean is true only if the double is != 0.
       */
      ++_b;
      while(
        _b < static_cast<VertexDescriptor>(_basePtr->_distances.outerSize())
        && (_basePtr->_distances(internal(_i), _b) == 0.0)
      ) {
        ++_b;
      }
      return *this;
    }

    inline in_group_edge_iterator operator ++ (int) {
      auto copy = *this;
      ++(*this);
      return copy;
    }

    inline bool operator == (const in_group_edge_iterator& other) {
      return (
        _isLeft == other._isLeft
        && _i == other._i
        && _b == other._b
        && _basePtr == other._basePtr
      );
    }
    inline bool operator != (const in_group_edge_iterator& other) {
      return !(*this == other);
    }

    inline EdgeDescriptor operator *() const {
      return {_i, target()};
    }

    inline VertexDescriptor target() const {
      if(_isLeft) {
        return left(_b);
      }

      return right(_b);
    }

    inline double weight() const {
      return _basePtr->upperBound(internal(_i), _b);
    }
  };

  inline in_group_edge_iterator in_group_edges_begin(const VertexDescriptor i) const {
    return {*this, i};
  }

  inline in_group_edge_iterator in_group_edges_end(const VertexDescriptor i) const {
    return {*this, i, true};
  }
};

} // namespace DistanceGeometry

} // namespace molassembler

} // namespace Scine

#endif
