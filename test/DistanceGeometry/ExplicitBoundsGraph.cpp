/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/graph/graph_concepts.hpp"
#include "boost/test/unit_test.hpp"

#include "Molassembler/DistanceGeometry/ExplicitBoundsGraph.h"
#include "Molassembler/DistanceGeometry/SpatialModel.h"
#include "Molassembler/IO.h"
#include "ShortestPathsGraphTests.h"

#include <iostream>
#include <iomanip>

inline std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

BOOST_AUTO_TEST_CASE(ExplicitBoundsGraphStructure, *boost::unit_test::label("DG")) {
  using namespace Scine::Molassembler;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator("stereocenter_detection_molecules")
  ) {
    Molecule molecule = IO::read(
      currentFilePath.string()
    );

    using EG = DistanceGeometry::ExplicitBoundsGraph;

    DistanceGeometry::SpatialModel spatialModel {molecule, DistanceGeometry::Configuration {}};

    EG explicitGraph {
      molecule.graph().inner(),
      spatialModel.makePairwiseBounds()
    };

    auto spg = explicitGraph.graph();

    ShortestPathsGraphConcept::checkEdgeStructure(spg, EG::left, EG::right);

    auto edge_iterators = boost::edges(spg);
    auto iter = edge_iterators.first;
    auto end = edge_iterators.second;

    ShortestPathsGraphConcept::checkEdgeIteration(iter, end, spg, EG::isLeft, EG::sameSide);

    // Generated distances matrix must satisfy triangle inequalities
    auto distancesMatrixResult = explicitGraph.makeDistanceMatrix(randomnessEngine());
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

    unsigned triangleInequalityFailures = 0;
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
            ++triangleInequalityFailures;
            std::cout << "The triangle inequality is falsified along i = " << i
              << ", j = " << j << ", k = " << k << "!"
              << " d(i, k) = " << (d(i, k)) << " > " << (d(i, j) + d(j, k))
              << " = d(i, j) + d(j, k)"
              << nl;
          }
        }
      }
    }

#ifdef MOLASSEMBLER_EXPLICIT_GRAPH_USE_SPECIALIZED_GOR1_ALGORITHM
    BOOST_CHECK_MESSAGE(
      static_cast<double>(triangleInequalityFailures) < std::max(1.0, 1e-3 * std::pow(matrN, 3)),
      "More than 0.1 % triangle inequality violations!"
        << nl << distancesMatrix << nl
    );
#else
    BOOST_CHECK_MESSAGE(
      triangleInequalityFailures == 0,
      "Generated distance matrix does not satisfy triangle inequalities!"
        << nl << distancesMatrix << nl
    );
  }
#endif
}
