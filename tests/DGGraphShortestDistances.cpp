#define BOOST_TEST_MODULE DGExplicitGraphTests
#define BOOST_TEST_DYN_LINK
#include "boost/test/unit_test.hpp"
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "DistanceGeometry/ImplicitGraphBoost.h"

#include "DistanceGeometry/MoleculeSpatialModel.h"
#include "DistanceGeometry/ExplicitGraph.h"
#include "DistanceGeometry/DistanceBoundsMatrix.h"
#include "template_magic/Numeric.h"
#include "template_magic/Enumerate.h"
#include "constexpr_magic/FloatingPointComparison.h"
#include "IO.h"

#include "boost/graph/bellman_ford_shortest_paths.hpp"
#include "Graph/Gor1.h"
#include "DistanceGeometry/Gor1.h"


#include <chrono>

/* TODO
 * VERY WEIRD: Algorithms executed on SPG find better shortest paths than on LG
 * or DBM (-> tighter bounds).
 */

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

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

// TODO Remove this, make a specialization of the Gor1Functor
std::vector<double> Gor1SPG (const MoleculeManip::DistanceGeometry::ImplicitGraph& graph, unsigned sourceVertex) {
  /* Prep */
  using Graph = MoleculeManip::DistanceGeometry::ImplicitGraph;
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
  const MoleculeManip::DistanceGeometry::DistanceBoundsMatrix& boundsRef;

  DBM_FW_Functor(const MoleculeManip::DistanceGeometry::DistanceBoundsMatrix& bounds) : boundsRef(bounds) {}

  MoleculeManip::DistanceGeometry::DistanceBoundsMatrix operator() () {
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

  using namespace MoleculeManip;

  boost::filesystem::path filesPath("../tests/mol_files/stereocenter_detection_molecules");
  boost::filesystem::recursive_directory_iterator end;

  IO::MOLFileHandler molHandler;
  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    Molecule molecule = molHandler.readSingle(
      currentFilePath.string()
    );

    DistanceGeometry::MoleculeSpatialModel spatialModel {
      molecule,
      DistanceGeometry::MoleculeSpatialModel::DistanceMethod::UFFLike
    };

    using LG = DistanceGeometry::ExplicitGraph;
    using Vertex = LG::GraphType::vertex_descriptor;

    LG lg {
      molecule,
      spatialModel.makeBoundList()
    };

    auto lgGraph = lg.getGraph();

    unsigned N = boost::num_vertices(lgGraph);

    for(Vertex a = 0; a < N / 2; ++a) {
      for(Vertex b = 0; b < N / 2; ++b) {
        if(a == b) {
          // No lara edge
          BOOST_CHECK_MESSAGE(
            !boost::edge(left(a), right(a), lgGraph).second,
            "Same-index l(a) -> r(a) exists for a = " << a
          );
        } else {
          auto lalb = boost::edge(left(a), left(b), lgGraph);
          if(lalb.second) {
            /* If there is a la -> lb edge, there must also be
             * - lb -> la (same weight)
             * - la -> rb
             * - lb -> ra
             * - ra -> rb (same weight)
             * - rb -> ra (same weight)
             */

            auto lalbWeight = boost::get(boost::edge_weight, lgGraph, lalb.first);

            BOOST_CHECK_MESSAGE(
              TemplateMagic::all_of(
                std::vector<decltype(lalb)> {
                  boost::edge(left(b), left(a), lgGraph),
                  boost::edge(right(b), right(a), lgGraph),
                  boost::edge(right(a), right(b), lgGraph)
                },
                [&](const auto& edgeFoundPair) -> bool {
                  if(!edgeFoundPair.second) {
                    // Edge must be found
                    return false;
                  }

                  auto weight = boost::get(boost::edge_weight, lgGraph, edgeFoundPair.first);
                  return weight == lalbWeight;
                }
              ),
              "Not all matching in-group edges were found or have identical weights for "
                << "a = " << a << ", b = " << b
            );

            BOOST_CHECK_MESSAGE(
              TemplateMagic::all_of(
                std::vector<decltype(lalb)> {
                  boost::edge(left(a), right(b), lgGraph),
                  boost::edge(left(b), right(a), lgGraph),
                },
                [&](const auto& edgeFoundPair) -> bool {
                  if(!edgeFoundPair.second) {
                    // Edge must be found
                    return false;
                  }

                  auto weight = boost::get(boost::edge_weight, lgGraph, edgeFoundPair.first);
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
  using namespace MoleculeManip;

  IO::MOLFileHandler molHandler;
  boost::filesystem::path filesPath("../tests/mol_files/stereocenter_detection_molecules");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;
    Molecule sampleMol = molHandler.readSingle(
      currentFilePath.string()
    );

    DistanceGeometry::MoleculeSpatialModel spatialModel {
      sampleMol,
      DistanceGeometry::MoleculeSpatialModel::DistanceMethod::UFFLike
    };

    const auto boundsList = spatialModel.makeBoundList();

    DistanceGeometry::ExplicitGraph limits {sampleMol, boundsList};
    DistanceGeometry::DistanceBoundsMatrix spatialModelBounds {sampleMol, boundsList};

    // This conforms to the triangle inequality bounds
    auto boundsMatrix = DBM_FW_Functor {spatialModelBounds} ();

    bool pass = true;
    for(unsigned a = 0; a < sampleMol.numAtoms(); ++a) {
      auto BF_LG_distances = BFFunctor<DistanceGeometry::ExplicitGraph::GraphType> {limits.getGraph()} (2 * a);

      auto Gor_LG_distances = Gor1Functor<DistanceGeometry::ExplicitGraph::GraphType> {limits.getGraph()} (2 * a);

      if(
        !TemplateMagic::all_of(
          TemplateMagic::zipMap(
            BF_LG_distances,
            Gor_LG_distances,
            [](const double& a, const double& b) -> bool {
              return ConstexprMagic::floating::isCloseRelative(a, b, 1e-8);
            }
          )
        )
      ) {
        pass = false;
        std::cout << "Not all pairs of Bellmann-Ford and Gor1 shortest-paths-distances on the ExplicitGraph are within 1e-8 relative tolerance!" << nl;
        break;
      }

      for(unsigned j = 0; j < BF_LG_distances.size(); j += 2) {
        if(j / 2 == a) {
          continue;
        }

        if(-BF_LG_distances.at(j + 1) > BF_LG_distances.at(j)) {
          pass = false;
          std::cout << "A lower bound is greater than the corresponding upper bound!" << nl;
          break;
        }
      }

      if(
        !TemplateMagic::all_of(
          enumerate(BF_LG_distances),
          [&boundsMatrix, &a](const auto& enumPair) -> bool {
            const auto& index = enumPair.index;
            const auto& distance = enumPair.value;

            if(index / 2 == a) {
              return true;
            }

            if(index % 2 == 0) {
              // To left index -> upper bound
              return ConstexprMagic::floating::isCloseRelative(
                boundsMatrix.upperBound(a, index / 2),
                distance,
                1e-4
              );
            }

            // To right index -> lower bound
            return ConstexprMagic::floating::isCloseRelative(
              boundsMatrix.lowerBound(a, index / 2),
              -distance,
              1e-4
            );
          }
        )
      ) {
        pass = false;
        std::cout << "Bellman-Ford ExplicitGraph shortest paths do not represent triangle inequality bounds!" << nl;
        std::cout << "Failed on a = " << a << nl;
        std::cout << "Distances:" << nl << TemplateMagic::condenseIterable(BF_LG_distances) << nl << boundsMatrix.access() << nl;
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
    
    using SPGVertex = DistanceGeometry::ImplicitGraph::VertexDescriptor;
    SPGVertex N = boost::num_vertices(shortestPathsGraph);

    // Check that both graphs are 1:1 identical
    bool identical = true;
    for(SPGVertex i = 0; i < N && identical; ++i) {
      for(SPGVertex j = 0; j < N; ++j) {
        auto spg_edge = boost::edge(i, j, shortestPathsGraph);
        auto lg_edge = boost::edge(i, j, limits.getGraph());

        if(spg_edge.second != lg_edge.second) {
          identical = false;
          std::cout << "Graphs do match for edge " << i << " -> " << j << ". "
            << "LG: " << lg_edge.second << ", SPG: " << spg_edge.second << nl;
        }

        if(spg_edge.second && lg_edge.second) {
          auto spg_edge_weight = boost::get(boost::edge_weight, shortestPathsGraph, spg_edge.first);
          auto lg_edge_weight = boost::get(boost::edge_weight, limits.getGraph(), lg_edge.first);

          if(spg_edge_weight != lg_edge_weight) {
            identical = false;
            std::cout << "Edge weight on edge " << i << " -> " << j << " does not match! "
              << "LG: " << lg_edge_weight << ", SPG: " << spg_edge_weight << nl;
          }
        }
      }
    }

    BOOST_REQUIRE_MESSAGE(identical, "LG and SPG are not identical graphs");

    pass = true;
    for(unsigned a = 0; a < sampleMol.numAtoms(); ++a) {
      // ImplicitGraph without implicit bounds should be consistent with ExplicitGraph
      auto BF_SPG_distances = BFFunctor<DistanceGeometry::ImplicitGraph> {shortestPathsGraph} (2 * a);

      auto Gor_SPG_distances = Gor1Functor<DistanceGeometry::ImplicitGraph> {shortestPathsGraph} (2 * a);

      if(
        !TemplateMagic::all_of(
          TemplateMagic::zipMap(
            BF_SPG_distances,
            Gor_SPG_distances,
            [](const double& a, const double& b) -> bool {
              return ConstexprMagic::floating::isCloseRelative(a, b, 1e-8);
            }
          )
        )
      ) {
        pass = false;
        std::cout << "Not all pairs of Bellmann-Ford and Gor1 shortest-paths-distances on the ImplicitGraph are within 1e-8 relative tolerance!" << nl;
        break;
      }
    
      auto spec_Gor_SPG_distances = Gor1SPG(shortestPathsGraph, 2 * a);

      if(
        !TemplateMagic::all_of(
          TemplateMagic::zipMap(
            Gor_SPG_distances,
            spec_Gor_SPG_distances,
            [](const double& a, const double& b) -> bool {
              return ConstexprMagic::floating::isCloseRelative(a, b, 1e-8);
            }
          )
        )
      ) {
        pass = false;
        std::cout << "Not all pairs of specialized and unspecialized Gor1 shortest-paths-distances on the ImplicitGraph are within 1e-8 relative tolerance!" << nl;
        std::cout << TemplateMagic::condenseIterable(spec_Gor_SPG_distances) << nl << nl
          << TemplateMagic::condenseIterable(Gor_SPG_distances) << nl << nl;
        break;
      }

      for(unsigned j = 0; j < Gor_SPG_distances.size(); j += 2) {
        if(j / 2 == a) {
          continue;
        }

        if(-Gor_SPG_distances.at(j + 1) > Gor_SPG_distances.at(j)) {
          pass = false;
          std::cout << "A lower bound is greater than the corresponding upper bound!" << nl;
          break;
        }
      }

      if(
        !TemplateMagic::all_of(
          enumerate(BF_SPG_distances),
          [&boundsMatrix, &a](const auto& enumPair) -> bool {
            const auto& index = enumPair.index;
            const auto& distance = enumPair.value;

            if(index / 2 == a) {
              return true;
            }

            if(index % 2 == 0) {
              // To left index -> upper bound
              return (
                ConstexprMagic::floating::isCloseRelative(
                  boundsMatrix.upperBound(a, index / 2),
                  distance,
                  1e-4
                ) || (
                  // Lower upper bounds are better
                  boundsMatrix.upperBound(a, index / 2) > distance
                )
              );
            }

            // To right index -> lower bound
            return ConstexprMagic::floating::isCloseRelative(
              boundsMatrix.lowerBound(a, index / 2),
              -distance,
              1e-4
            ) || (
              // Larger lower bounds are better
              boundsMatrix.lowerBound(a, index / 2) < -distance
            );
          }
        )
      ) {
        pass = false;
        std::cout << "Bellman-Ford ImplicitGraph shortest paths do not represent triangle inequality bounds!" << nl;
        std::cout << "Failed on a = " << a << ", j = {";
        for(unsigned j = 0; j < BF_SPG_distances.size(); ++j) {
          if(j / 2 != a) {
            if(j % 2 == 0) {
              if(
                !ConstexprMagic::floating::isCloseRelative(
                  boundsMatrix.upperBound(a, j / 2),
                  BF_SPG_distances.at(j),
                  1e-4
                )
              ) {
                std::cout << j << ", ";
              }
            } else {
              if(
                !ConstexprMagic::floating::isCloseRelative(
                  boundsMatrix.lowerBound(a, j / 2),
                  -BF_SPG_distances.at(j),
                  1e-4
                )
              ) {
                std::cout << j << ", ";
              }
            }
          }
        }
        std::cout << "}" << nl;
        std::cout << "Distances:" << nl << TemplateMagic::condenseIterable(BF_SPG_distances) << nl << boundsMatrix.access() << nl;
        break;
      }
    }

    BOOST_CHECK_MESSAGE(
      pass,
      "Bellman-Ford and Gor ImplicitGraph shortest paths distances fail consistency checks with Floyd-Warshall DistanceBoundsMatrix."
    );
  }
}
