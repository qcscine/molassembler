/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_TEST_MODULE DGStrainedMolecules
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "Molassembler/Conformers.h"
#include "Molassembler/DistanceGeometry/SpatialModel.h"
#include "Molassembler/IO.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/StereopermutatorList.h"

#include <iostream>

using namespace Scine;
using namespace Molassembler;

BOOST_AUTO_TEST_CASE(transSpanningImpossibilitiesRemoved) {
  auto mol = IO::read("inorganics/multidentate/Co(ox)3.mol");

  const auto stereopermutatorOption = mol.stereopermutators().option(0);
  BOOST_REQUIRE(stereopermutatorOption);

  const unsigned N = stereopermutatorOption -> numAssignments();

  for(unsigned i = 0; i < N; ++i) {
    mol.assignStereopermutator(0, i);

    auto ensemble = generateRandomEnsemble(mol, 10);
    for(auto& positionResult : ensemble) {
      if(!positionResult) {
        BOOST_FAIL(positionResult.error().message());
      }
    }
  }
}

void readFileGenConformationAndWriteFile(const boost::filesystem::path& filePath) {
  using namespace std::string_literals;

  std::cout << "Processing " << filePath.stem().string() << std::endl;

  // Read the file
  auto mol = IO::read(filePath.string());

  DistanceGeometry::SpatialModel spatialModel {mol, DistanceGeometry::Configuration {}};

  spatialModel.writeGraphviz(filePath.stem().string() + ".dot"s);

  try {
    // Generate a conformation
    if(auto positionsResult = generateRandomConformation(mol)) {
      // Write the generated conformation to file
      IO::write(
        filePath.stem().string() + "-generated.mol"s,
        mol,
        positionsResult.value()
      );
    } else {
      BOOST_FAIL(positionsResult.error().message());
    }
  } catch(std::exception& e) {
    std::cout << "Unhandled exception: " << e.what() << std::endl;
    BOOST_FAIL("Unhandled exception encountered in conformation generation");
    throw;
  }
}

BOOST_AUTO_TEST_CASE(strainedOrganicMolecules) {
  boost::filesystem::path filesPath("strained_organic_molecules");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;
    BOOST_CHECK_NO_THROW(
      readFileGenConformationAndWriteFile(currentFilePath)
    );
  }
}
