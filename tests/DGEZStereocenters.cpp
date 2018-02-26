#define BOOST_TEST_MODULE DGEZStereocenters
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include <iostream>
#include "IO.h"
#include "DistanceGeometry/generateConformation.h"

void readFileGenConformationAndWriteFile(const boost::filesystem::path& filePath) {
  using namespace MoleculeManip;

  std::cout << "Processing " << filePath.stem().string() << std::endl;

  // Read the file
  IO::MOLFileHandler molHandler;
  auto mol = molHandler.readSingle(filePath.string());

  std::cout << mol << std::endl;

  // Generate a conformation
  auto positionsResult = DistanceGeometry::generateConformation(mol);
  BOOST_REQUIRE_MESSAGE(positionsResult, positionsResult.error().message());
  auto positions = positionsResult.value();

  // Write the generated conformation to file
  molHandler.writeSingle(
    filePath.stem().string() + "-generated.mol"s,
    mol,
    positions
  );
}

BOOST_AUTO_TEST_CASE(strainedOrganicMolecules) {
  boost::filesystem::path filesPath("../tests/mol_files/ez_stereocenters");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;
    BOOST_CHECK_NO_THROW(
      readFileGenConformationAndWriteFile(currentFilePath)
    );
  }
}
