#define BOOST_TEST_MODULE ImplicitGraphTestModule
#define BOOST_TEST_DYN_LINK
#include "boost/test/unit_test.hpp"
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "DistanceGeometry/ImplicitGraphBoost.h"
#include "DistanceGeometry/MoleculeSpatialModel.h"

#include "boost/graph/graph_concepts.hpp"
#include "IO.h"

#include <iostream>
#include <iomanip>

inline std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

BOOST_AUTO_TEST_CASE(graphConcepts) {
  using namespace molassembler;

  using GraphType = molassembler::DistanceGeometry::ImplicitGraph;

  BOOST_CONCEPT_ASSERT(( boost::VertexListGraphConcept<GraphType> ));
  BOOST_CONCEPT_ASSERT(( boost::EdgeListGraphConcept<GraphType> ));
  BOOST_CONCEPT_ASSERT(( boost::AdjacencyMatrixConcept<GraphType> ));

  BOOST_CONCEPT_ASSERT(( 
    boost::ReadablePropertyMapConcept<
      boost::property_map<GraphType, boost::vertex_index_t>::type,
      boost::graph_traits<GraphType>::vertex_descriptor
    >
  ));

  BOOST_CONCEPT_ASSERT(( 
    boost::ReadablePropertyMapConcept<
      boost::property_map<GraphType, boost::edge_weight_t>::type,
      boost::graph_traits<GraphType>::edge_descriptor
    >
  ));

  /*BOOST_CONCEPT_ASSERT(( 
    boost::WritablePropertyMapConcept<
      boost::property_map<GraphType, boost::edge_weight_t>::type,
      boost::graph_traits<GraphType>::edge_descriptor
    >
  ));*/

  /*BOOST_CONCEPT_ASSERT(( 
    boost::PropertyGraphConcept<
      GraphType,
      boost::graph_traits<GraphType>::edge_descriptor,
      boost::edge_weight_t
    >
  ));*/

  BOOST_CONCEPT_ASSERT(( boost::IncidenceGraphConcept<GraphType> ));
}

template<typename UnsignedType>
UnsignedType left(UnsignedType a) {
  return 2 * a;
}

template<typename UnsignedType>
UnsignedType right(UnsignedType a) {
  return 2 * a + 1;
}

BOOST_AUTO_TEST_CASE(nonVisualTests) {
  using namespace molassembler;

  boost::filesystem::path filesPath("../tests/mol_files/stereocenter_detection_molecules");
  boost::filesystem::recursive_directory_iterator end;

  IO::MOLFileHandler molHandler;
  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;

    Molecule molecule = molHandler.readSingle(
      currentFilePath.string()
    );

    using SPG = DistanceGeometry::ImplicitGraph;

    DistanceGeometry::MoleculeSpatialModel spatialModel {
      molecule,
      DistanceGeometry::MoleculeSpatialModel::DistanceMethod::UFFLike
    };

    SPG spg {
      molecule,
      spatialModel.makeBoundList()
    };

    SPG::VertexDescriptor N = boost::num_vertices(spg);

    for(SPG::VertexDescriptor i = 0; i < N; ++i) {
      BOOST_CHECK_MESSAGE(
        static_cast<SPG::VertexDescriptor>(
          std::distance(
            spg.obegin(i),
            spg.oend(i)
          )
        ) == boost::out_degree(i, spg),
        "Out degree of vertex " << i << " does not match out_edge_iterator begin-end distance"
      );

      // Out-edge iterator tests
      auto iter = spg.obegin(i);
      auto end = spg.oend(i);

      while(iter != end) {
        auto edgeDescriptor = *iter;

        BOOST_CHECK_MESSAGE(
          boost::target(edgeDescriptor, spg) == iter.target(),
          "out_edge_iter boost-target and iter-member-target do not match for {"
          << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
        );

        BOOST_CHECK_MESSAGE(
          boost::get(boost::edge_weight, spg, edgeDescriptor) == iter.weight(),
          "out_edge_iter boost-get-weight and iter-member-weight do not match for {"
          << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
        );

        ++iter;
      }

      // in_group_edge_iterator tests
      auto groupIter = spg.in_group_edges_begin(i);
      const auto groupEnd = spg.in_group_edges_end(i);

      while(groupIter != groupEnd) {
        auto edgeDescriptor = *groupIter;
        auto target = groupIter.target();
        auto weight = groupIter.weight();
        auto boostWeight = boost::get(boost::edge_weight, spg, edgeDescriptor);

        BOOST_CHECK_MESSAGE(
          boost::edge(i, target, spg).second,
          "Edge reported by in_group_edge_iterator does not exist per boost::edge. "
            << "iterator: " << i << " -> " << target
        );

        BOOST_CHECK_MESSAGE(
          weight == boostWeight,
          "in-group-edge does not yield same weight via boost-get and weight method"
        );

        BOOST_CHECK_MESSAGE(
          weight > 0,
          "in-group-edge weight isn't greater than zero for edge " << i << " -> " 
            << target << ", weight = " << weight
        );

        ++groupIter;
      }
    }

    for(SPG::VertexDescriptor a = 0; a < N / 2; ++a) {
      BOOST_CHECK_MESSAGE(
        !boost::edge(2 * a, 2 * a + 1, spg).second,
        "Same-a edge exists for a = " << a
      );

      for(SPG::VertexDescriptor b = 0; b < N / 2; ++b) {
        if(a == b) {
          continue;
        }

        // No right-to-left edges
        BOOST_CHECK_MESSAGE(
          !boost::edge(right(a), left(b), spg).second,
          "r(a) -> l(b) for a = " << a << ", b = " << b
        );

        BOOST_CHECK_MESSAGE(
          !boost::edge(right(b), left(a), spg).second,
          "r(b) -> l(a) for a = " << a << ", b = " << b
        );

        // If there is an edge from la to lb
        auto lalb = boost::edge(2 * a, 2 * b, spg);
        if(lalb.second) {
          auto lalbWeight = boost::get(boost::edge_weight, spg, lalb.first);

          auto lbra = boost::edge(2 * b, 2 * a + 1, spg);

          BOOST_CHECK_MESSAGE(
            lbra.second,
            "matching l(b) -> r(a) edge does not exist for a = " << a
              << ", b = " << b
          );

          auto lbraWeight = boost::get(boost::edge_weight, spg, lbra.first);

          auto larb = boost::edge(2 * a, 2 * b + 1, spg);

          BOOST_CHECK_MESSAGE(
            larb.second,
            "matching l(a) -> r(b) edge does not exist for a = " << a
              << ", b = " << b
          );

          auto larbWeight = boost::get(boost::edge_weight, spg, larb.first);

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

    BOOST_CHECK_MESSAGE(
      static_cast<SPG::VertexDescriptor>(
        std::distance(
          spg.ebegin(),
          spg.eend()
        ) 
      ) == boost::num_edges(spg),
      "Number of edges does not match edge iterator begin-end distance"
    );

    auto iter = spg.ebegin();
    auto end = spg.eend();

    while(iter != end) {
      auto edgeDescriptor = *iter;

      auto boostSource = boost::source(edgeDescriptor, spg);
      auto boostTarget = boost::target(edgeDescriptor, spg);

      // Boost and iter target fetch match
      BOOST_CHECK_MESSAGE(
        boostTarget == iter.target(),
        "edge_iter boost-target and iter-member-target do not match for {"
        << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
      );

      // No right-to-left edges
      BOOST_CHECK_MESSAGE(
        !(!SPG::isLeft(boostSource) && SPG::isLeft(boostTarget)),
        "Edge points from right to left! " << boostSource << " -> " << boostTarget
      );

      auto boostEdgeWeight = boost::get(boost::edge_weight, spg, edgeDescriptor);

      // Boost and iter weight fetch match
      BOOST_CHECK_MESSAGE(
        boostEdgeWeight == iter.weight(),
        "edge_iter boost-get-weight and iter-member-weight do not match for {"
        << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
      );

      if(boostSource % 2 == boostTarget % 2) {
        // In-group edge

        // Reverse exists
        auto reverseEdge = boost::edge(boostTarget, boostSource, spg);
        BOOST_CHECK_MESSAGE(
          reverseEdge.second, 
          "Reverse edge does not exist for in-group edge " 
            << boostSource << " -> " << boostTarget
        );

        // Reverse has same edge weight
        BOOST_CHECK_MESSAGE(
          boost::get(boost::edge_weight, spg, reverseEdge.first) == boostEdgeWeight,
          "Reverse edge for " << boostSource << " -> " << boostTarget 
            << " does not have same edge weight"
        );
      }

      ++iter;
    }

    // Generated distances matrix must satisfy triangle inequalities
    auto distancesMatrixResult = spg.makeDistanceMatrix();
    BOOST_REQUIRE_MESSAGE(distancesMatrixResult, distancesMatrixResult.error().message());
    auto distancesMatrix = distancesMatrixResult.value();

    auto d = [&distancesMatrix](const AtomIndexType& i, const AtomIndexType& j) -> double {
      return distancesMatrix(
        std::min(i, j),
        std::max(i, j)
      );
    };

    bool passTriangleInequalities = true;
    unsigned matrN = distancesMatrix.cols();
    for(AtomIndexType i = 0; i < matrN; ++i) {
      for(AtomIndexType j = 0; j < matrN; ++j) {
        if(i == j) {
          continue;
        }

        for(AtomIndexType k = 0; k < matrN; ++k) {
          if(j == k) {
            continue;
          }

          if(d(i, k) > d(i, j) + d(j, k)) {
            passTriangleInequalities = false;
            std::cout << "The triangle inequality is falsified along i = " << i
              << ", j = " << j << ", k = " << k << "!" 
              << " d(i, k) = " << (d(i, k)) << " > " << (d(i, j) + d(j, k))
              << " = d(i, j) + d(j, k)"
              << nl;
          }
        }
      }
    }

    BOOST_CHECK_MESSAGE(
      passTriangleInequalities,
      "Generated distance matrix does not satisfy triangle inequalities!"
        << nl << distancesMatrix << nl
    );
  }
}
