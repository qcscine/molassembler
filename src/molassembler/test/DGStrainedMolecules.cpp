#define BOOST_TEST_MODULE DGStrainedMolecules
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include <iostream>
#include "IO.h"
#include "DistanceGeometry/generateConformation.h"

BOOST_AUTO_TEST_CASE(transSpanningImpossibilitiesRemoved) {
  using namespace molassembler;

  auto mol = IO::read("test_files/inorganics/multidentate/Co(ox)3.mol");

  const auto& stereocenterPtr = mol.getStereocenterList().at(0);
  unsigned N = stereocenterPtr -> numStereopermutations();

  for(unsigned i = 0; i < N; ++i) {
    mol.assignStereocenter(0, i);

    auto ensembleResult = DistanceGeometry::generateEnsemble(mol, 10);
    if(!ensembleResult) {
      BOOST_FAIL(ensembleResult.error().message());
    }
  }
}

void readFileGenConformationAndWriteFile(const boost::filesystem::path& filePath) {
  using namespace molassembler;

  std::cout << "Processing " << filePath.stem().string() << std::endl;

  // Read the file
  auto mol = IO::read(filePath.string());

  DistanceGeometry::MoleculeSpatialModel spatialModel {mol};

  spatialModel.writeGraphviz(filePath.stem().string() + ".dot"s);

  try {
    // Generate a conformation
    if(auto positionsResult = DistanceGeometry::generateConformation(mol)) {
      // Write the generated conformation to file
      IO::write(
        filePath.stem().string() + "-generated.mol"s,
        mol,
        positionsResult.value()
      );
    } else {
      BOOST_FAIL(positionsResult.error().message());
    }
  } catch(std::exception e) {
    std::cout << "Unhandled exception: " << e.what() << std::endl;
    BOOST_FAIL("Unhandled exception encountered in conformation generation");
    throw;
  }
}

BOOST_AUTO_TEST_CASE(strainedOrganicMolecules) {
  boost::filesystem::path filesPath("test_files/strained_organic_molecules");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;
    BOOST_CHECK_NO_THROW(
      readFileGenConformationAndWriteFile(currentFilePath)
    );
  }
}
