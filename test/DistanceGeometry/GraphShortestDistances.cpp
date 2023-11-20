/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/graph/two_bit_color_map.hpp"
#include "boost/test/unit_test.hpp"

#include "Molassembler/Temple/Adaptors/Enumerate.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/constexpr/FloatingPointComparison.h"
#include "Molassembler/Temple/constexpr/Numeric.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"

// DO NOT CHANGE THIS INCLUDE ORDER (implicit graph needs to go first)
#include "Molassembler/DistanceGeometry/ImplicitBoundsGraphBoost.h"

#include "Molassembler/DistanceGeometry/DistanceBoundsMatrix.h"
#include "Molassembler/DistanceGeometry/ExplicitBoundsGraph.h"
#include "Molassembler/DistanceGeometry/SpatialModel.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/IO.h"
#include "Molassembler/Molecule.h"

#include "Molassembler/Graph/Gor1.h"
#include "Molassembler/DistanceGeometry/Gor1.h"


// This include order may seem weird, but it is necessary like this
#include "boost/graph/bellman_ford_shortest_paths.hpp"

#include <chrono>

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

using namespace Scine::Molassembler;
using namespace DistanceGeometry;

template<class Graph>
struct BfFunctor {
  const Graph& graph;

  BfFunctor(const Graph& g) : graph(g) {}

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

// Use specialized GOR1 implementation for ImplicitBoundsGraph
std::vector<double> Gor1IG (const DistanceGeometry::ImplicitBoundsGraph& graph, unsigned sourceVertex) {
  /* Prep */
  using Graph = DistanceGeometry::ImplicitBoundsGraph;
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
  const DistanceGeometry::DistanceBoundsMatrix& boundsRef;

  DBM_FW_Functor(const DistanceGeometry::DistanceBoundsMatrix& bounds) : boundsRef(bounds) {}

  DistanceGeometry::DistanceBoundsMatrix operator() () {
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

BOOST_AUTO_TEST_CASE(ShortestPathsGraphConcepts, *boost::unit_test::label("DG")) {

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("stereocenter_detection_molecules")
  ) {
    Molecule molecule = IO::read(
      currentFilePath.string()
    );

    DistanceGeometry::SpatialModel spatialModel {molecule, DistanceGeometry::Configuration {}};

    using EG = DistanceGeometry::ExplicitBoundsGraph;
    using Vertex = EG::GraphType::vertex_descriptor;

    EG eg {
      molecule.graph().inner(),
      spatialModel.makePairwiseBounds()
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
              Temple::all_of(
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
              Temple::all_of(
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

bool shortestPathsRepresentInequalities(
  const ExplicitBoundsGraph& limits,
  const DistanceBoundsMatrix& boundsMatrix,
  const unsigned N
) {
  bool pass = true;
  for(unsigned outerVertex = 0; outerVertex < N; ++outerVertex) {
    auto BF_EG_distances = BfFunctor<ExplicitBoundsGraph::GraphType> {limits.graph()} (2 * outerVertex);
    auto Gor_EG_distances = Gor1Functor<ExplicitBoundsGraph::GraphType> {limits.graph()} (2 * outerVertex);

    if(
      !Temple::all_of(
        Temple::Adaptors::zip(
          BF_EG_distances,
          Gor_EG_distances
        ),
        [](const double a, const double b) -> bool {
          return Temple::Floating::isCloseRelative(a, b, 1e-8);
        }
      )
    ) {
      pass = false;
      std::cout << "Not all pairs of Bellmann-Ford and Gor1 shortest-paths-distances on the ExplicitBoundsGraph are within 1e-8 relative tolerance!" << nl;
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
      !Temple::all_of(
        Temple::Adaptors::enumerate(BF_EG_distances),
        [&boundsMatrix, &outerVertex](const auto& enumPair) -> bool {
          const auto& index = enumPair.index;
          const auto& distance = enumPair.value;

          if(index / 2 == outerVertex) {
            return true;
          }

          if(index % 2 == 0) {
            // To left index -> upper bound
            return Temple::Floating::isCloseRelative(
              boundsMatrix.upperBound(outerVertex, index / 2),
              distance,
              1e-4
            );
          }

          // To right index -> lower bound
          return Temple::Floating::isCloseRelative(
            boundsMatrix.lowerBound(outerVertex, index / 2),
            -distance,
            1e-4
          );
        }
      )
    ) {
      pass = false;
      std::cout << "Bellman-Ford ExplicitBoundsGraph shortest paths do not represent triangle inequality bounds!" << nl;
      std::cout << "Failed on outerVertex = " << outerVertex << nl;
      std::cout << "Distances:" << nl << Temple::condense(BF_EG_distances) << nl << boundsMatrix.access() << nl;
      break;
    }
  }
  return pass;
}

bool graphsIdentical(const ExplicitBoundsGraph& explicitGraph, const ImplicitBoundsGraph& implicitGraph) {
  using IgVertex = ImplicitBoundsGraph::VertexDescriptor;
  IgVertex N = boost::num_vertices(implicitGraph);

  // Check that both graphs are 1:1 identical
  bool identical = true;
  for(IgVertex i = 0; i < N && identical; ++i) {
    for(IgVertex j = 0; j < N; ++j) {
      auto igEdge = boost::edge(i, j, implicitGraph);
      auto egEdge = boost::edge(i, j, explicitGraph.graph());

      if(igEdge.second != egEdge.second) {
        identical = false;
        std::cout << "Graphs do not match for edge " << i << " -> " << j << ". "
          << "EG: " << egEdge.second << ", IG: " << igEdge.second << nl;
      }

      if(igEdge.second && egEdge.second) {
        auto ig_edge_weight = boost::get(boost::edge_weight, implicitGraph, igEdge.first);
        auto eg_edge_weight = boost::get(boost::edge_weight, explicitGraph.graph(), egEdge.first);

        if(ig_edge_weight != eg_edge_weight) {
          identical = false;
          std::cout << "Edge weight on edge " << i << " -> " << j << " does not match! "
            << "EG: " << eg_edge_weight << ", IG: " << ig_edge_weight << nl;
        }
      }
    }
  }

  return identical;
}

bool shortestGraphsAlgorithmsResultsMatch(
  const unsigned N,
  const ImplicitBoundsGraph& implicitGraph,
  const DistanceBoundsMatrix& boundsMatrix
) {
  bool pass = true;
  for(unsigned outerVertex = 0; outerVertex < N; ++outerVertex) {
    // ImplicitBoundsGraph without implicit bounds should be consistent with ExplicitBoundsGraph
    auto BF_IG_distances = BfFunctor<DistanceGeometry::ImplicitBoundsGraph> {implicitGraph} (2 * outerVertex);

    auto Gor_IG_distances = Gor1Functor<DistanceGeometry::ImplicitBoundsGraph> {implicitGraph} (2 * outerVertex);

    if(
      !Temple::all_of(
        Temple::Adaptors::zip(BF_IG_distances, Gor_IG_distances),
        [](const double a, const double b) -> bool {
          return Temple::Floating::isCloseRelative(a, b, 1e-8);
        }
      )
    ) {
      pass = false;
      std::cout << "Not all pairs of Bellmann-Ford and Gor1 shortest-paths-distances on the ImplicitBoundsGraph are within 1e-8 relative tolerance!" << nl;
      break;
    }

    auto spec_Gor_IG_distances = Gor1IG(implicitGraph, 2 * outerVertex);

    if(
      !Temple::all_of(
        Temple::Adaptors::zip(Gor_IG_distances, spec_Gor_IG_distances),
        [](const double a, const double b) -> bool {
          return Temple::Floating::isCloseRelative(a, b, 1e-8);
        }
      )
    ) {
      pass = false;
      std::cout << "Not all pairs of specialized and unspecialized Gor1 shortest-paths-distances on the ImplicitBoundsGraph are within 1e-8 relative tolerance!" << nl;
      std::cout << Temple::condense(spec_Gor_IG_distances) << nl << nl
        << Temple::condense(Gor_IG_distances) << nl << nl;
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
      !Temple::all_of(
        Temple::Adaptors::enumerate(BF_IG_distances),
        [&boundsMatrix, &outerVertex](const auto& enumPair) -> bool {
          const auto& index = enumPair.index;
          const auto& distance = enumPair.value;

          if(index / 2 == outerVertex) {
            return true;
          }

          if(index % 2 == 0) {
            // To left index -> upper bound
            return (
              Temple::Floating::isCloseRelative(
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
          return Temple::Floating::isCloseRelative(
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
      std::cout << "Bellman-Ford ImplicitBoundsGraph shortest paths do not represent triangle inequality bounds!" << nl;
      std::cout << "Failed on outerVertex = " << outerVertex << ", j = {";
      for(unsigned j = 0; j < BF_IG_distances.size(); ++j) {
        if(j / 2 != outerVertex) {
          if(j % 2 == 0) {
            if(
              !Temple::Floating::isCloseRelative(
                boundsMatrix.upperBound(outerVertex, j / 2),
                BF_IG_distances.at(j),
                1e-4
              )
            ) {
              std::cout << j << ", ";
            }
          } else {
            if(
              !Temple::Floating::isCloseRelative(
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
      std::cout << "Distances:" << nl << Temple::condense(BF_IG_distances) << nl << boundsMatrix.access() << nl;
      break;
    }
  }

  return pass;
}

BOOST_AUTO_TEST_CASE(GraphShortestDistancesCorrectness, *boost::unit_test::label("DG")) {
  /* Want to do combination tests with the iodo-alkane.
   *
   * Can calculate shortest paths with either:
   * - Floyd-Warshall in DistanceBoundsMatrix (all-pairs!)
   * - Bellman-Ford
   * - Simplified Gor1
   *
   * On either:
   * - ExplicitBoundsGraph (Underlying fully explicit BGL graph)
   * - ImplicitBoundsGraph (underlying ValueBounds, but BGL interface)
   *
   * Time and compare correctness:
   * - DistanceBoundsMatrix  + Floyd-Warshall (currently ONLY correct impl)
   * - ExplicitBoundsGraph + Bellman-Ford
   * - ExplicitBoundsGraph + Gor1
   * - ImplicitBoundsGraph + Bellman-Ford
   * - ImplicitBoundsGraph + Gor1
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
  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("stereocenter_detection_molecules")
  ) {
    Molecule sampleMol = IO::read(
      currentFilePath.string()
    );
    const unsigned N = sampleMol.graph().V();

    DistanceGeometry::SpatialModel spatialModel {sampleMol, DistanceGeometry::Configuration {}};

    const auto boundsList = spatialModel.makePairwiseBounds();

    DistanceGeometry::ExplicitBoundsGraph explicitGraph {sampleMol.graph().inner(), boundsList};
    DistanceGeometry::DistanceBoundsMatrix spatialModelBounds {sampleMol.graph().inner(), boundsList};

    // This conforms to the triangle inequality bounds
    auto boundsMatrix = DBM_FW_Functor {spatialModelBounds} ();

    BOOST_CHECK_MESSAGE(
      shortestPathsRepresentInequalities(explicitGraph, boundsMatrix, N),
      "Bellman-Ford and/or Gor1 ExplicitBoundsGraph shortest-paths fails consistency tests with Floyd-Warshall DistanceBoundsMatrix."
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
    DistanceGeometry::ImplicitBoundsGraph implicitGraph {sampleMol.graph().inner(), boundsList};

    BOOST_REQUIRE_MESSAGE(
      graphsIdentical(explicitGraph, implicitGraph),
      "EG and IG are not identical graphs"
    );

    BOOST_CHECK_MESSAGE(
      shortestGraphsAlgorithmsResultsMatch(N, implicitGraph, boundsMatrix),
      "Bellman-Ford and Gor ImplicitBoundsGraph shortest paths distances fail consistency checks with Floyd-Warshall DistanceBoundsMatrix."
    );
  }
}
