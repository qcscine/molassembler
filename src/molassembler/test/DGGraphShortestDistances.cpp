/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/graph/two_bit_color_map.hpp"
#include "boost/test/unit_test.hpp"

#include "temple/Adaptors/Enumerate.h"
#include "temple/Adaptors/Zip.h"
#include "temple/Functional.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Stringify.h"

#include "molassembler/DistanceGeometry/ImplicitGraphBoost.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/DistanceGeometry/ExplicitGraph.h"
#include "molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule.h"

#include "gor1/Gor1.h"
#include "molassembler/DistanceGeometry/Gor1.h"

// This include order may seem weird, but it is necessary like this
#include "boost/graph/bellman_ford_shortest_paths.hpp"

#include <chrono>

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

using namespace Scine;

template<class Graph>
struct BFFunctor {
  const Graph& graph;

  BFFunctor(const Graph& g) : graph(g) {}

  std::vector<double> operator() (unsigned sourceVertex) {
    unsigned M = boost::num_vertices(graph);

    std::vector<double> distance (M);

    boost::bellman_ford_shortest_paths(
      graph,
      M,
      boost::root_vertex(sourceVertex).
      distance_map(&distance[0])
    );

    return distance;
  }
};

template<class Graph>
struct Gor1Functor {
  const Graph& graph;

  Gor1Functor(const Graph& g) : graph(g) {}

  std::vector<double> operator() (unsigned sourceVertex) {
    /* Prep */
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;

    unsigned N = boost::num_vertices(graph);
    std::vector<double> distances(N);
    std::vector<Vertex> predecessors(N);

    auto predecessor_map = boost::make_iterator_property_map(
      predecessors.begin(),
      get(boost::vertex_index, graph)
    );

    auto distance_map = boost::make_iterator_property_map(
      distances.begin(),
      get(boost::vertex_index, graph)
    );

    using ColorMapType = boost::two_bit_color_map<>;
    ColorMapType color_map {N};

    /* Execution */

    boost::gor1_simplified_shortest_paths(
      graph,
      Vertex {sourceVertex},
      predecessor_map,
      color_map,
      distance_map
    );

    return distances;
  }
};

// Use specialized GOR1 implementation for ImplicitGraph
std::vector<double> Gor1IG (const molassembler::DistanceGeometry::ImplicitGraph& graph, unsigned sourceVertex) {
  /* Prep */
  using Graph = molassembler::DistanceGeometry::ImplicitGraph;
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;

  unsigned N = boost::num_vertices(graph);
  std::vector<double> distances(N);
  std::vector<Vertex> predecessors(N);

  auto predecessor_map = boost::make_iterator_property_map(
    predecessors.begin(),
    get(boost::vertex_index, graph)
  );

  auto distance_map = boost::make_iterator_property_map(
    distances.begin(),
    get(boost::vertex_index, graph)
  );

  using ColorMapType = boost::two_bit_color_map<>;
  ColorMapType color_map {N};

  /* Execution */

  boost::gor1_ig_shortest_paths(
    graph,
    Vertex {sourceVertex},
    predecessor_map,
    color_map,
    distance_map
  );

  return distances;
}

struct DBM_FW_Functor {
  const molassembler::DistanceGeometry::DistanceBoundsMatrix& boundsRef;

  DBM_FW_Functor(const molassembler::DistanceGeometry::DistanceBoundsMatrix& bounds) : boundsRef(bounds) {}

  molassembler::DistanceGeometry::DistanceBoundsMatrix operator() () {
    auto boundsCopy = boundsRef;

    boundsCopy.smooth();

    return boundsCopy;
  }
};

template<typename UnsignedType>
UnsignedType left(UnsignedType a) {
  return 2 * a;
}

template<typename UnsignedType>
UnsignedType right(UnsignedType a) {
  return 2 * a + 1;
}

BOOST_AUTO_TEST_CASE(conceptTests) {
  using namespace Scine;
  using namespace molassembler;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("stereocenter_detection_molecules")
  ) {
    Molecule molecule = IO::read(
      currentFilePath.string()
    );

    DistanceGeometry::SpatialModel spatialModel {molecule, DistanceGeometry::Configuration {}};

    using EG = DistanceGeometry::ExplicitGraph;
    using Vertex = EG::GraphType::vertex_descriptor;

    EG eg {
      molecule,
      spatialModel.makeBoundsList()
    };

    auto egGraph = eg.graph();

    unsigned N = boost::num_vertices(egGraph);

    for(Vertex a = 0; a < N / 2; ++a) {
      for(Vertex b = 0; b < N / 2; ++b) {
        if(a == b) {
          // No lara edge
          BOOST_CHECK_MESSAGE(
            !boost::edge(left(a), right(a), egGraph).second,
            "Same-index l(a) -> r(a) exists for a = " << a
          );
        } else {
          auto lalb = boost::edge(left(a), left(b), egGraph);
          if(lalb.second) {
            /* If there is a la -> lb edge, there must also be
             * - lb -> la (same weight)
             * - la -> rb
             * - lb -> ra
             * - ra -> rb (same weight)
             * - rb -> ra (same weight)
             */

            auto lalbWeight = boost::get(boost::edge_weight, egGraph, lalb.first);

            BOOST_CHECK_MESSAGE(
              temple::all_of(
                std::vector<decltype(lalb)> {
                  boost::edge(left(b), left(a), egGraph),
                  boost::edge(right(b), right(a), egGraph),
                  boost::edge(right(a), right(b), egGraph)
                },
                [&](const auto& edgeFoundPair) -> bool {
                  if(!edgeFoundPair.second) {
                    // Edge must be found
                    return false;
                  }

                  auto weight = boost::get(boost::edge_weight, egGraph, edgeFoundPair.first);
                  return weight == lalbWeight;
                }
              ),
              "Not all matching in-group edges were found or have identical weights for "
                << "a = " << a << ", b = " << b
            );

            BOOST_CHECK_MESSAGE(
              temple::all_of(
                std::vector<decltype(lalb)> {
                  boost::edge(left(a), right(b), egGraph),
                  boost::edge(left(b), right(a), egGraph),
                },
                [&](const auto& edgeFoundPair) -> bool {
                  if(!edgeFoundPair.second) {
                    // Edge must be found
                    return false;
                  }

                  auto weight = boost::get(boost::edge_weight, egGraph, edgeFoundPair.first);
                  return std::fabs(weight) < lalbWeight;
                }
              ),
              "Not all matching cross-group edges were found or have lower absolute weights for "
                << "a = " << a << ", b = " << b
            );
          }
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(correctnessTests) {

  /* Want to do combination tests with the iodo-alkane.
   *
   * Can calculate shortest paths with either:
   * - Floyd-Warshall in DistanceBoundsMatrix (all-pairs!)
   * - Bellman-Ford
   * - Simplified Gor1
   *
   * On either:
   * - ExplicitGraph (Underlying fully explicit BGL graph)
   * - ImplicitGraph (underlying ValueBounds, but BGL interface)
   *
   * Time and compare correctness:
   * - DistanceBoundsMatrix  + Floyd-Warshall (currently ONLY correct impl)
   * - ExplicitGraph + Bellman-Ford
   * - ExplicitGraph + Gor1
   * - ImplicitGraph + Bellman-Ford
   * - ImplicitGraph + Gor1
   *
   * The bottom four combinations should yield equal shortest paths distances,
   * although they are not consistent with the triangle inequalities. One such
   * shortest paths calculation yields only the triangle inequality limits for
   * one pair of atom indices.
   *
   * Only DBM + FW is assumed correct.
   */
  using namespace std::chrono;

  // DBM + FW
  using namespace molassembler;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("stereocenter_detection_molecules")
  ) {
    Molecule sampleMol = IO::read(
      currentFilePath.string()
    );

    DistanceGeometry::SpatialModel spatialModel {sampleMol, DistanceGeometry::Configuration {}};

    const auto boundsList = spatialModel.makeBoundsList();

    DistanceGeometry::ExplicitGraph limits {sampleMol, boundsList};
    DistanceGeometry::DistanceBoundsMatrix spatialModelBounds {sampleMol, boundsList};

    // This conforms to the triangle inequality bounds
    auto boundsMatrix = DBM_FW_Functor {spatialModelBounds} ();

    bool pass = true;
    for(unsigned outerVertex = 0; outerVertex < sampleMol.graph().N(); ++outerVertex) {
      auto BF_EG_distances = BFFunctor<DistanceGeometry::ExplicitGraph::GraphType> {limits.graph()} (2 * outerVertex);

      auto Gor_EG_distances = Gor1Functor<DistanceGeometry::ExplicitGraph::GraphType> {limits.graph()} (2 * outerVertex);

      if(
        !temple::all_of(
          temple::adaptors::zip(
            BF_EG_distances,
            Gor_EG_distances
          ),
          [](const double a, const double b) -> bool {
            return temple::floating::isCloseRelative(a, b, 1e-8);
          }
        )
      ) {
        pass = false;
        std::cout << "Not all pairs of Bellmann-Ford and Gor1 shortest-paths-distances on the ExplicitGraph are within 1e-8 relative tolerance!" << nl;
        break;
      }

      for(unsigned j = 0; j < BF_EG_distances.size(); j += 2) {
        if(j / 2 == outerVertex) {
          continue;
        }

        if(-BF_EG_distances.at(j + 1) > BF_EG_distances.at(j)) {
          pass = false;
          std::cout << "A lower bound is greater than the corresponding upper bound!" << nl;
          break;
        }
      }

      if(
        !temple::all_of(
          temple::adaptors::enumerate(BF_EG_distances),
          [&boundsMatrix, &outerVertex](const auto& enumPair) -> bool {
            const auto& index = enumPair.index;
            const auto& distance = enumPair.value;

            if(index / 2 == outerVertex) {
              return true;
            }

            if(index % 2 == 0) {
              // To left index -> upper bound
              return temple::floating::isCloseRelative(
                boundsMatrix.upperBound(outerVertex, index / 2),
                distance,
                1e-4
              );
            }

            // To right index -> lower bound
            return temple::floating::isCloseRelative(
              boundsMatrix.lowerBound(outerVertex, index / 2),
              -distance,
              1e-4
            );
          }
        )
      ) {
        pass = false;
        std::cout << "Bellman-Ford ExplicitGraph shortest paths do not represent triangle inequality bounds!" << nl;
        std::cout << "Failed on outerVertex = " << outerVertex << nl;
        std::cout << "Distances:" << nl << temple::condense(BF_EG_distances) << nl << boundsMatrix.access() << nl;
        break;
      }
    }

    BOOST_CHECK_MESSAGE(
      pass,
      "Bellman-Ford and/or Gor1 ExplicitGraph shortest-paths fails consistency tests with Floyd-Warshall DistanceBoundsMatrix."
    );

    /* If the shortest paths calculation can yield the full upper and lower
     * bounds within the triangle inequality limits for one source atom, then the
     * full triangle inequality bounds can be calculated within
     * O(N * O(shortest paths calc)). While generating the distances matrix,
     * this full calculation is unnecessary, a single shortest paths calculation
     * yields the bounds for the random distance choice. The full complexity of
     * creating a distance matrix should be O(N² * O(shortest paths calc)).
     *
     * Hopefully, the complexity of the shortest paths calculation can stay
     * constant as distance bounds are added / narrowed. This is difficult, and
     * perhaps direct access to the emerging distance matrix is necessary to
     * ensure O(1) bounds access!
     */
    DistanceGeometry::ImplicitGraph shortestPathsGraph {sampleMol, boundsList};

    using IGVertex = DistanceGeometry::ImplicitGraph::VertexDescriptor;
    IGVertex N = boost::num_vertices(shortestPathsGraph);

    // Check that both graphs are 1:1 identical
    bool identical = true;
    for(IGVertex i = 0; i < N && identical; ++i) {
      for(IGVertex j = 0; j < N; ++j) {
        auto ig_edge = boost::edge(i, j, shortestPathsGraph);
        auto eg_edge = boost::edge(i, j, limits.graph());

        if(ig_edge.second != eg_edge.second) {
          identical = false;
          std::cout << "Graphs do not match for edge " << i << " -> " << j << ". "
            << "EG: " << eg_edge.second << ", IG: " << ig_edge.second << nl;
        }

        if(ig_edge.second && eg_edge.second) {
          auto ig_edge_weight = boost::get(boost::edge_weight, shortestPathsGraph, ig_edge.first);
          auto eg_edge_weight = boost::get(boost::edge_weight, limits.graph(), eg_edge.first);

          if(ig_edge_weight != eg_edge_weight) {
            identical = false;
            std::cout << "Edge weight on edge " << i << " -> " << j << " does not match! "
              << "EG: " << eg_edge_weight << ", IG: " << ig_edge_weight << nl;
          }
        }
      }
    }

    BOOST_REQUIRE_MESSAGE(identical, "EG and IG are not identical graphs");

    pass = true;
    for(unsigned outerVertex = 0; outerVertex < sampleMol.graph().N(); ++outerVertex) {
      // ImplicitGraph without implicit bounds should be consistent with ExplicitGraph
      auto BF_IG_distances = BFFunctor<DistanceGeometry::ImplicitGraph> {shortestPathsGraph} (2 * outerVertex);

      auto Gor_IG_distances = Gor1Functor<DistanceGeometry::ImplicitGraph> {shortestPathsGraph} (2 * outerVertex);

      if(
        !temple::all_of(
          temple::adaptors::zip(
            BF_IG_distances,
            Gor_IG_distances
          ),
          [](const double a, const double b) -> bool {
            return temple::floating::isCloseRelative(a, b, 1e-8);
          }
        )
      ) {
        pass = false;
        std::cout << "Not all pairs of Bellmann-Ford and Gor1 shortest-paths-distances on the ImplicitGraph are within 1e-8 relative tolerance!" << nl;
        break;
      }

      auto spec_Gor_IG_distances = Gor1IG(shortestPathsGraph, 2 * outerVertex);

      if(
        !temple::all_of(
          temple::adaptors::zip(
            Gor_IG_distances,
            spec_Gor_IG_distances
          ),
          [](const double a, const double b) -> bool {
            return temple::floating::isCloseRelative(a, b, 1e-8);
          }
        )
      ) {
        pass = false;
        std::cout << "Not all pairs of specialized and unspecialized Gor1 shortest-paths-distances on the ImplicitGraph are within 1e-8 relative tolerance!" << nl;
        std::cout << temple::condense(spec_Gor_IG_distances) << nl << nl
          << temple::condense(Gor_IG_distances) << nl << nl;
        break;
      }

      for(unsigned j = 0; j < Gor_IG_distances.size(); j += 2) {
        if(j / 2 == outerVertex) {
          continue;
        }

        if(-Gor_IG_distances.at(j + 1) > Gor_IG_distances.at(j)) {
          pass = false;
          std::cout << "A lower bound is greater than the corresponding upper bound!" << nl;
          break;
        }
      }

      if(
        !temple::all_of(
          temple::adaptors::enumerate(BF_IG_distances),
          [&boundsMatrix, &outerVertex](const auto& enumPair) -> bool {
            const auto& index = enumPair.index;
            const auto& distance = enumPair.value;

            if(index / 2 == outerVertex) {
              return true;
            }

            if(index % 2 == 0) {
              // To left index -> upper bound
              return (
                temple::floating::isCloseRelative(
                  boundsMatrix.upperBound(outerVertex, index / 2),
                  distance,
                  1e-4
                ) || (
                  // Lower upper bounds are better
                  boundsMatrix.upperBound(outerVertex, index / 2) > distance
                )
              );
            }

            // To right index -> lower bound
            return temple::floating::isCloseRelative(
              boundsMatrix.lowerBound(outerVertex, index / 2),
              -distance,
              1e-4
            ) || (
              // Larger lower bounds are better
              boundsMatrix.lowerBound(outerVertex, index / 2) < -distance
            );
          }
        )
      ) {
        pass = false;
        std::cout << "Bellman-Ford ImplicitGraph shortest paths do not represent triangle inequality bounds!" << nl;
        std::cout << "Failed on outerVertex = " << outerVertex << ", j = {";
        for(unsigned j = 0; j < BF_IG_distances.size(); ++j) {
          if(j / 2 != outerVertex) {
            if(j % 2 == 0) {
              if(
                !temple::floating::isCloseRelative(
                  boundsMatrix.upperBound(outerVertex, j / 2),
                  BF_IG_distances.at(j),
                  1e-4
                )
              ) {
                std::cout << j << ", ";
              }
            } else {
              if(
                !temple::floating::isCloseRelative(
                  boundsMatrix.lowerBound(outerVertex, j / 2),
                  -BF_IG_distances.at(j),
                  1e-4
                )
              ) {
                std::cout << j << ", ";
              }
            }
          }
        }
        std::cout << "}" << nl;
        std::cout << "Distances:" << nl << temple::condense(BF_IG_distances) << nl << boundsMatrix.access() << nl;
        break;
      }
    }

    BOOST_CHECK_MESSAGE(
      pass,
      "Bellman-Ford and Gor ImplicitGraph shortest paths distances fail consistency checks with Floyd-Warshall DistanceBoundsMatrix."
    );
  }
}
