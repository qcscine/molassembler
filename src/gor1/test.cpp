/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_TEST_MODULE Gor1TestsModule
#include "boost/test/unit_test.hpp"
#include "boost/graph/random.hpp"
#include "boost/graph/bellman_ford_shortest_paths.hpp"
#include "boost/graph/two_bit_color_map.hpp"

#include "gor1/Gor1.h"

#include <random>

using EdgeWeightProperty = boost::property<boost::edge_weight_t, double>;

using Graph = boost::adjacency_list<
  boost::setS,
  boost::vecS,
  boost::directedS,
  boost::no_property,
  EdgeWeightProperty
>;

class RNG {
private:
  std::vector<unsigned> _seeds;

  void _initGenerator() {

#ifdef NDEBUG
    std::random_device randomDevice;
    for(unsigned n = 0; n < 5; n++) _seeds.emplace_back(randomDevice());
#else
    _seeds.emplace_back(2721813754);
#endif

    std::seed_seq _seedSequence(_seeds.begin(), _seeds.end());
    engine.seed(_seedSequence);
  }

public:
  mutable std::mt19937 engine;

  RNG() {
    _initGenerator();
  }

  double get(double lower, double upper) const {
    assert(lower <= upper);
    std::uniform_real_distribution<double> uniformDistribution(lower, upper);
    return uniformDistribution(engine);
  }
};

static RNG rng;

// Imitate the testing setup in the paper (see the Gor1 header)
Graph SPACYC_P2N(
  const unsigned vertices,
  const unsigned edges,
  const double lower,
  const double upper
) {
  Graph graph {vertices};

  // Set up the main path
  for(unsigned i = 0; i < vertices - 1; ++i) {
    auto edgeAddPair = boost::add_edge(i, i + 1, graph);
    assert(edgeAddPair.second);
    boost::get(boost::edge_weight, graph, edgeAddPair.first) = rng.get(lower, upper);
  }

  // Add random edges
  for(unsigned N = 0; N < edges; ++N) {
    Graph::vertex_descriptor i, j;

    // Keep fetching new pairs of indices if i is j or i and j are already connected
    do {
      i = boost::random_vertex(graph, rng.engine);
      j = boost::random_vertex(graph, rng.engine);
    } while(
      i == j
      || boost::edge(std::min(i, j), std::max(i, j), graph).second
    );

    auto edgeAddPair = boost::add_edge(std::min(i, j), std::max(i, j), graph);
    boost::get(boost::edge_weight, graph, edgeAddPair.first) = rng.get(lower, upper);
  }

  return graph;
}

BOOST_AUTO_TEST_CASE(gor1Tests) {
  unsigned V = 8192;
  unsigned E = 131072;
  double l = -10000;
  double u = 10000;

  for(unsigned i = 0; i < 10; ++i) {
    Graph testGraph = SPACYC_P2N(V, E, l, u);

    /* GOR1 */
    std::vector<double> distances (V);
    std::vector<Graph::vertex_descriptor> predecessors (V);
    using ColorMapType = boost::two_bit_color_map<>;
    ColorMapType color_map {V};

    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      boost::get(boost::vertex_index, testGraph)
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      boost::get(boost::vertex_index, testGraph)
    );

    // fill color map with white
    std::fill(
      color_map.data.get(),
      color_map.data.get() + (color_map.n + ColorMapType::elements_per_char - 1)
        / ColorMapType::elements_per_char,
      0
    );

    boost::gor1_simplified_shortest_paths(
      testGraph,
      0,
      predecessor_map,
      color_map,
      distance_map
    );

    /* Bellman-Ford */
    std::vector<double> BFdistances (V);
    std::vector<Graph::vertex_descriptor> BFpredecessors (V);

    auto BFpredecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      boost::get(boost::vertex_index, testGraph)
    );

    auto BFdistance_map = boost::make_iterator_property_map(
      BFdistances.begin(),
      boost::get(boost::vertex_index, testGraph)
    );

    boost::bellman_ford_shortest_paths(
      testGraph,
      V,
      boost::predecessor_map(BFpredecessor_map).
      distance_map(BFdistance_map).
      root_vertex(0)
    );

    BOOST_CHECK_MESSAGE(
      distances == BFdistances,
      "Calculated shortest paths are not the same between BF and GOR1"
    );
  }
}
