#define BOOST_TEST_MODULE RingDecomposition
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "IO.h"

void decomposeMol(const MoleculeManip::Molecule& mol) {
  auto cycleData = mol.getAdjacencyList().getCycleData();

  for(
    auto cycleIter = cycleData.getIterator();
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    const auto cycleEdges = cycleIter.getCurrentCycle();
    std::cout << cycleEdges.size() << std::endl;
  }
}

void readAndDecompose(const boost::filesystem::path& filePath) {
  using namespace MoleculeManip;

  std::cout << "Processing " << filePath.stem().string() << std::endl;

  // Read the file
  IO::MOLFileHandler molHandler;
  auto mol = molHandler.readSingle(filePath.string());

  decomposeMol(mol);
}

BOOST_AUTO_TEST_CASE(ringDecomposition) {
  boost::filesystem::path filesPath("../tests/mol_files/strained_organic_molecules");
  boost::filesystem::recursive_directory_iterator end;

  for(boost::filesystem::recursive_directory_iterator i(filesPath); i != end; i++) {
    const boost::filesystem::path currentFilePath = *i;
    BOOST_CHECK_NO_THROW(
      readAndDecompose(currentFilePath)
    );
  }
}
