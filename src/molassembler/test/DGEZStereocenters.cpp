#define BOOST_TEST_MODULE DGBondStereocenters
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include <iostream>
#include "molassembler/IO.h"
#include "molassembler/Conformers.h"

void readFileGenConformationAndWriteFile(const boost::filesystem::path& filePath) {
  using namespace molassembler;

  std::cout << "Processing " << filePath.stem().string() << std::endl;

  // Read the file
  auto mol = IO::read(filePath.string());

  std::cout << mol << std::endl;

  // Generate a conformation
  auto positionsResult = generateConformation(mol);

  if(!positionsResult) {
    BOOST_FAIL(positionsResult.error().message());
  }

  auto positions = positionsResult.value();

  // Write the generated conformation to file
  IO::write(
    filePath.stem().string() + "-generated.mol"s,
    mol,
    positions
  );
}

BOOST_AUTO_TEST_CASE(ezStereocenterExamples) {
  boost::filesystem::path filesPath("test_files/ez_stereocenters");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;
    BOOST_CHECK_NO_THROW(
      readFileGenConformationAndWriteFile(currentFilePath)
    );
  }
}
