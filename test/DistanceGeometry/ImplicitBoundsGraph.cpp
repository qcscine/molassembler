/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"

#include "molassembler/DistanceGeometry/ImplicitBoundsGraphBoost.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"

#include "boost/graph/graph_concepts.hpp"
#include "molassembler/IO.h"
#include "ShortestPathsGraphTests.h"

#include <iostream>
#include <iomanip>

inline std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

BOOST_AUTO_TEST_CASE(ImplicitBoundsGraphConcepts) {
  using namespace Scine::Molassembler;

  using GraphType = DistanceGeometry::ImplicitBoundsGraph;

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

BOOST_AUTO_TEST_CASE(ImplicitBoundsGraphStructure) {
  using namespace Scine::Molassembler;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("stereocenter_detection_molecules")
  ) {
    Molecule molecule = IO::read(
      currentFilePath.string()
    );

    using IG = DistanceGeometry::ImplicitBoundsGraph;

    DistanceGeometry::SpatialModel spatialModel {molecule, DistanceGeometry::Configuration {}};

    IG ig {
      molecule.graph().inner(),
      spatialModel.makePairwiseBounds()
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

    ShortestPathsGraphConcept::checkEdgeStructure(ig, IG::left, IG::right);

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

    ShortestPathsGraphConcept::checkEdgeIteration(iter, end, ig, IG::isLeft, IG::sameSide);

    // Generated distances matrix must satisfy triangle inequalities
    auto distancesMatrixResult = ig.makeDistanceMatrix(randomnessEngine());
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

#ifdef MOLASSEMBLER_IMPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
    BOOST_CHECK_MESSAGE(
      static_cast<double>(triangleInequalitiesFailures) < std::max(1.0, 1e-3 * std::pow(matrN, 3)),
      "More than 0.1 % triangle inequality violations!"
        << nl << distancesMatrix << nl
        << "In matrix of size " << matrN << ", got "
        << triangleInequalitiesFailures
        << " triangle inequality violations. There are "
        << std::pow(matrN, 3) << " such relations in this matrix. ("
        << (100 * static_cast<double>(triangleInequalitiesFailures) / std::pow(matrN, 3))
        << " %)" << nl
    );
#else
    BOOST_CHECK_MESSAGE(
      triangleInequalitiesFailures == 0,
      "There are triangle inequality violations!" << nl << distancesMatrix << nl
    );
#endif
  }
}
