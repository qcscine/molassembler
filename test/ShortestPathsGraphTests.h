/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#ifndef INCLUDE_MOLASSEMBLER_TESTS_SHORTEST_PATHS_GRAPH_TESTS_H
#define INCLUDE_MOLASSEMBLER_TESTS_SHORTEST_PATHS_GRAPH_TESTS_H

/* Helper file to test shortest paths graph concepts
 */

namespace detail {

template<class>
struct sfinae_true : std::true_type {};

template<class Iterator>
static auto testHasTarget(int) -> sfinae_true<
  decltype(std::declval<Iterator>().target())
>;

template<class Iterator>
static auto testHasTarget(long) -> std::false_type;

template<class Iterator>
struct hasTarget : decltype(testHasTarget<Iterator>(0)) {};

} // namespace detail

struct ShortestPathsGraphConcept {
  template<typename Graph, typename LeftFunctor, typename RightFunctor>
  static void checkEdgeStructure(
    const Graph& g,
    LeftFunctor&& left,
    RightFunctor&& right
  ) {
    /* For every atomic index a ∈ [0, N), the graph contains two vertices,
     * left(a) and right(a). For every pair of atomic indices a, b that have a
     * distance bound, i.e. a lower and upper bound on their spatial distance
     * within the current metric space, there are six directed edges in the
     * graph:
     *
     * - Two edges (one bidirectional edge) between left(a) and left(b) with the
     *   edge weight of the upper bound
     * - Two edges between right(a) and right(b) with the edge weight of
     *   the upper bound
     * - A unidirectional edge from left(a) to right(b) with the edge weight of
     *   the lower bound, negated
     * - A unidirectional edge from left(b) to right(a) with the edge weight of
     *   the lower bound, negated
     */

    using VertexDescriptor = decltype(boost::num_vertices(std::declval<Graph>()));
    const VertexDescriptor N = boost::num_vertices(g);

    for(VertexDescriptor a = 0; a < N / 2; ++a) {
      BOOST_CHECK_MESSAGE(
        !boost::edge(left(a), right(a), g).second,
        "Same-a edge exists for a = " << a
      );

      for(VertexDescriptor b = 0; b < N / 2; ++b) {
        if(a == b) {
          continue;
        }

        // No right-to-left edges
        BOOST_CHECK_MESSAGE(
          !boost::edge(right(a), left(b), g).second,
          "r(a) -> l(b) for a = " << a << ", b = " << b
        );

        BOOST_CHECK_MESSAGE(
          !boost::edge(right(b), left(a), g).second,
          "r(b) -> l(a) for a = " << a << ", b = " << b
        );

        // If there is an edge from la to lb
        auto lalb = boost::edge(left(a), left(b), g);
        if(lalb.second) {
          auto lalbWeight = boost::get(boost::edge_weight, g, lalb.first);

          auto lbra = boost::edge(left(b), right(a), g);

          BOOST_CHECK_MESSAGE(
            lbra.second,
            "matching l(b) -> r(a) edge does not exist for a = " << a
              << ", b = " << b
          );

          auto lbraWeight = boost::get(boost::edge_weight, g, lbra.first);

          auto larb = boost::edge(left(a), right(b), g);

          BOOST_CHECK_MESSAGE(
            larb.second,
            "matching l(a) -> r(b) edge does not exist for a = " << a
              << ", b = " << b
          );

          auto larbWeight = boost::get(boost::edge_weight, g, larb.first);

          BOOST_CHECK_MESSAGE(
            larbWeight == lbraWeight,
            "l(a) -> r(b) and l(b) -> r(a) edges do not have same weight for a = "
              << a << ", b = " << b << ": " << larbWeight << ", " << lbraWeight
          );

          BOOST_CHECK_MESSAGE(
            lalbWeight > std::fabs(larbWeight),
            "l(a) -> l(b) weight isn't greater than abs of l(a) -> r(b) weight for a = "
              << a << ", b = " << b
          );
        }
      }
    }
  }

  template<typename Iter, typename VertexDescriptor>
  static std::enable_if_t<detail::hasTarget<Iter>::value, void> checkIterShortcuts(
    const Iter& iter,
    const VertexDescriptor source,
    const VertexDescriptor target,
    const double weight
  ) {
    // Boost and iter target fetch match
    BOOST_CHECK_MESSAGE(
      target == iter.target(),
      "edge_iter boost-target and iter-member-target do not match for {"
      << source << ", " << target << " != " << iter.target() << "}"
    );

    BOOST_CHECK_MESSAGE(
      weight == iter.weight(),
      "edge_iter boost-get-weight and iter-member-weight do not match for {"
      << source << ", " << iter.target() << "}"
    );
  }

  template<typename Iter, typename VertexDescriptor>
  static std::enable_if_t<!detail::hasTarget<Iter>::value, void> checkIterShortcuts(
    const Iter& /* iter */,
    const VertexDescriptor /* source */,
    const VertexDescriptor /* target */,
    const double /* weight */
  ) {}

  template<typename Iter, typename Graph, typename IsLeftFunctor, typename SameSideFunctor>
  static void checkEdgeIteration(
    Iter iter,
    Iter end,
    const Graph& g,
    IsLeftFunctor&& isLeft,
    SameSideFunctor&& sameSide
  ) {
    while(iter != end) {
      auto edgeDescriptor = *iter;

      auto boostSource = boost::source(edgeDescriptor, g);
      auto boostTarget = boost::target(edgeDescriptor, g);

      // No right-to-left edges
      BOOST_CHECK_MESSAGE(
        !(!isLeft(boostSource) && isLeft(boostTarget)),
        "Edge points from right to left! " << boostSource << " -> " << boostTarget
      );

      auto boostEdgeWeight = boost::get(boost::edge_weight, g, edgeDescriptor);

      checkIterShortcuts(iter, boostSource, boostTarget, boostEdgeWeight);

      if(sameSide(boostSource, boostTarget)) {
        // Reverse exists
        auto reverseEdge = boost::edge(boostTarget, boostSource, g);
        BOOST_CHECK_MESSAGE(
          reverseEdge.second,
          "Reverse edge does not exist for in-group edge "
            << boostSource << " -> " << boostTarget
        );

        // Reverse has same edge weight
        BOOST_CHECK_MESSAGE(
          boost::get(boost::edge_weight, g, reverseEdge.first) == boostEdgeWeight,
          "Reverse edge for " << boostSource << " -> " << boostTarget
            << " does not have same edge weight"
        );
      }

      ++iter;
    }
  }
};

#endif
