/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_TEST_MODULE DGStrainedMolecules
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "molassembler/Conformers.h"
#include "molassembler/DistanceGeometry/SpatialModel.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule.h"
#include "molassembler/StereopermutatorList.h"

#include <iostream>

using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(transSpanningImpossibilitiesRemoved) {
  auto mol = IO::read("inorganics/multidentate/Co(ox)3.mol");

  const auto stereopermutatorOption = mol.stereopermutators().option(0);
  BOOST_REQUIRE(stereopermutatorOption);

  const unsigned N = stereopermutatorOption -> numAssignments();

  for(unsigned i = 0; i < N; ++i) {
    mol.assignStereopermutator(0, i);

    auto ensemble = generateEnsemble(mol, 10);
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

  DistanceGeometry::SpatialModel spatialModel {mol, DistanceGeometry::Configuration {}, randomnessEngine()};

  spatialModel.writeGraphviz(filePath.stem().string() + ".dot"s);

  try {
    // Generate a conformation
    if(auto positionsResult = generateConformation(mol)) {
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
