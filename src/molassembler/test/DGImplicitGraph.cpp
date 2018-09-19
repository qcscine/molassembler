// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#define BOOST_TEST_MODULE ImplicitGraphTestModule
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"

#include "molassembler/DistanceGeometry/ImplicitGraphBoost.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"

#include "boost/graph/graph_concepts.hpp"
#include "molassembler/IO.h"

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

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("stereocenter_detection_molecules")
  ) {
    Molecule molecule = IO::read(
      currentFilePath.string()
    );

    using IG = DistanceGeometry::ImplicitGraph;

    DistanceGeometry::SpatialModel spatialModel {molecule};

    IG ig {
      molecule,
      spatialModel.makeBoundsList()
    };

    IG::VertexDescriptor N = boost::num_vertices(ig);

    for(IG::VertexDescriptor innerVertex = 0; innerVertex < N; ++innerVertex) {
      BOOST_CHECK_MESSAGE(
        static_cast<IG::VertexDescriptor>(
          std::distance(
            ig.obegin(innerVertex),
            ig.oend(innerVertex)
          )
        ) == boost::out_degree(innerVertex, ig),
        "Out degree of vertex " << innerVertex << " does not match out_edge_iterator begin-end distance"
      );

      // Out-edge iterator tests
      auto iter = ig.obegin(innerVertex);
      auto end = ig.oend(innerVertex);

      while(iter != end) {
        auto edgeDescriptor = *iter;

        BOOST_CHECK_MESSAGE(
          boost::target(edgeDescriptor, ig) == iter.target(),
          "out_edge_iter boost-target and iter-member-target do not match for {"
          << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
        );

        BOOST_CHECK_MESSAGE(
          boost::get(boost::edge_weight, ig, edgeDescriptor) == iter.weight(),
          "out_edge_iter boost-get-weight and iter-member-weight do not match for {"
          << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
        );

        ++iter;
      }

      // in_group_edge_iterator tests
      auto groupIter = ig.in_group_edges_begin(innerVertex);
      const auto groupEnd = ig.in_group_edges_end(innerVertex);

      while(groupIter != groupEnd) {
        auto edgeDescriptor = *groupIter;
        auto target = groupIter.target();
        auto weight = groupIter.weight();
        auto boostWeight = boost::get(boost::edge_weight, ig, edgeDescriptor);

        BOOST_CHECK_MESSAGE(
          boost::edge(innerVertex, target, ig).second,
          "Edge reported by in_group_edge_iterator does not exist per boost::edge. "
            << "iterator: " << innerVertex << " -> " << target
        );

        BOOST_CHECK_MESSAGE(
          weight == boostWeight,
          "in-group-edge does not yield same weight via boost-get and weight method"
        );

        BOOST_CHECK_MESSAGE(
          weight > 0,
          "in-group-edge weight isn't greater than zero for edge " << innerVertex << " -> "
            << target << ", weight = " << weight
        );

        ++groupIter;
      }
    }

    for(IG::VertexDescriptor a = 0; a < N / 2; ++a) {
      BOOST_CHECK_MESSAGE(
        !boost::edge(2 * a, 2 * a + 1, ig).second,
        "Same-a edge exists for a = " << a
      );

      for(IG::VertexDescriptor b = 0; b < N / 2; ++b) {
        if(a == b) {
          continue;
        }

        // No right-to-left edges
        BOOST_CHECK_MESSAGE(
          !boost::edge(right(a), left(b), ig).second,
          "r(a) -> l(b) for a = " << a << ", b = " << b
        );

        BOOST_CHECK_MESSAGE(
          !boost::edge(right(b), left(a), ig).second,
          "r(b) -> l(a) for a = " << a << ", b = " << b
        );

        // If there is an edge from la to lb
        auto lalb = boost::edge(2 * a, 2 * b, ig);
        if(lalb.second) {
          auto lalbWeight = boost::get(boost::edge_weight, ig, lalb.first);

          auto lbra = boost::edge(2 * b, 2 * a + 1, ig);

          BOOST_CHECK_MESSAGE(
            lbra.second,
            "matching l(b) -> r(a) edge does not exist for a = " << a
              << ", b = " << b
          );

          auto lbraWeight = boost::get(boost::edge_weight, ig, lbra.first);

          auto larb = boost::edge(2 * a, 2 * b + 1, ig);

          BOOST_CHECK_MESSAGE(
            larb.second,
            "matching l(a) -> r(b) edge does not exist for a = " << a
              << ", b = " << b
          );

          auto larbWeight = boost::get(boost::edge_weight, ig, larb.first);

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
      static_cast<IG::VertexDescriptor>(
        std::distance(
          ig.ebegin(),
          ig.eend()
        )
      ) == boost::num_edges(ig),
      "Number of edges does not match edge iterator begin-end distance"
    );

    auto iter = ig.ebegin();
    auto end = ig.eend();

    while(iter != end) {
      auto edgeDescriptor = *iter;

      auto boostSource = boost::source(edgeDescriptor, ig);
      auto boostTarget = boost::target(edgeDescriptor, ig);

      // Boost and iter target fetch match
      BOOST_CHECK_MESSAGE(
        boostTarget == iter.target(),
        "edge_iter boost-target and iter-member-target do not match for {"
        << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
      );

      // No right-to-left edges
      BOOST_CHECK_MESSAGE(
        !(!IG::isLeft(boostSource) && IG::isLeft(boostTarget)),
        "Edge points from right to left! " << boostSource << " -> " << boostTarget
      );

      auto boostEdgeWeight = boost::get(boost::edge_weight, ig, edgeDescriptor);

      // Boost and iter weight fetch match
      BOOST_CHECK_MESSAGE(
        boostEdgeWeight == iter.weight(),
        "edge_iter boost-get-weight and iter-member-weight do not match for {"
        << edgeDescriptor.first << ", " << edgeDescriptor.second << "}"
      );

      if(boostSource % 2 == boostTarget % 2) {
        // In-group edge

        // Reverse exists
        auto reverseEdge = boost::edge(boostTarget, boostSource, ig);
        BOOST_CHECK_MESSAGE(
          reverseEdge.second,
          "Reverse edge does not exist for in-group edge "
            << boostSource << " -> " << boostTarget
        );

        // Reverse has same edge weight
        BOOST_CHECK_MESSAGE(
          boost::get(boost::edge_weight, ig, reverseEdge.first) == boostEdgeWeight,
          "Reverse edge for " << boostSource << " -> " << boostTarget
            << " does not have same edge weight"
        );
      }

      ++iter;
    }

    // Generated distances matrix must satisfy triangle inequalities
    auto distancesMatrixResult = ig.makeDistanceMatrix();
    if(!distancesMatrixResult) {
      BOOST_FAIL(distancesMatrixResult.error().message());
    }
    auto distancesMatrix = distancesMatrixResult.value();

    auto d = [&distancesMatrix](const AtomIndex i, const AtomIndex j) -> double {
      return distancesMatrix(
        std::min(i, j),
        std::max(i, j)
      );
    };

    unsigned triangleInequalitiesFailures = 0;
    unsigned matrN = distancesMatrix.cols();
    for(AtomIndex i = 0; i < matrN; ++i) {
      for(AtomIndex j = 0; j < matrN; ++j) {
        if(i == j) {
          continue;
        }

        for(AtomIndex k = 0; k < matrN; ++k) {
          if(j == k) {
            continue;
          }

          if(d(i, k) > d(i, j) + d(j, k)) {
            ++triangleInequalitiesFailures;
            std::cout << "The triangle inequality is falsified along i = " << i
              << ", j = " << j << ", k = " << k << "!"
              << " d(i, k) = " << (d(i, k)) << " > " << (d(i, j) + d(j, k))
              << " = d(i, j) + d(j, k)"
              << nl;
          }
        }
      }
    }

    std::cout << "In matrix of size " << matrN << ", got "
      << triangleInequalitiesFailures
      << " triangle inequality violations. There are "
      << std::pow(matrN, 3) << " such relations in this matrix. ("
      << (100 * static_cast<double>(triangleInequalitiesFailures) / std::pow(matrN, 3))
      << " %)" << nl;

#ifdef MOLASSEMBLER_IMPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
    BOOST_CHECK_MESSAGE(
      static_cast<double>(triangleInequalitiesFailures) < std::max(1.0, 1e-3 * std::pow(matrN, 3)),
      "More than 0.1 % triangle inequality violations!"
        << nl << distancesMatrix << nl
    );
#else
    BOOST_CHECK_MESSAGE(
      triangleInequalitiesFailures == 0,
      "There are triangle inequality violations!" << nl << distancesMatrix << nl
    );
#endif
  }
}
