#define BOOST_TEST_MODULE RingDecomposition
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "template_magic/Containers.h"

#include "IO.h"

void testCycles(const MoleculeManip::CycleData& cycleData) {

  // Check that maxCycleSize works properly
  std::vector<unsigned> cycleSizes;

  for(
    auto cycleIter = cycleData.getCyclesIterator();
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    const auto cycleEdges = cycleIter.getCurrentCycle();
  // for(const auto& cycleEdges : cycleData) {
    cycleSizes.emplace_back(cycleEdges.size());
  }

  std::vector<unsigned> cyclesSmallerThanSix;

  for(
    auto cycleIter = cycleData.getCyclesIteratorSizeLE(5);
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    const auto cycleEdges = cycleIter.getCurrentCycle();
  // for(const auto& cycleEdges : cycleData.iterateCyclesSmallerThan(6)) {
    cyclesSmallerThanSix.emplace_back(cycleEdges.size());
  }

  BOOST_CHECK(!cycleSizes.empty());

  BOOST_CHECK(cycleSizes.size() >= cyclesSmallerThanSix.size());

  BOOST_CHECK(
    TemplateMagic::all_of(
      TemplateMagic::map(
        cyclesSmallerThanSix,
        [](const unsigned& cycleSize) -> bool {
          return cycleSize <= 6;
        }
      )
    )
  );

  std::cout << TemplateMagic::condenseIterable(cycleSizes) << std::endl;

  // Check that checking for cycles that contain a particular atom works
  /*for(
    auto cycleIter = cycleData.getCyclesIteratorContaining(5ul);
    !cycleIter.atEnd();
    cycleIter.advance()
  ) {
    const auto cycleEdges = cycleIter.getCurrentCycle();
    cycleSizes.emplace_back(cycleEdges.size());
  }*/
}

void readAndDecompose(const boost::filesystem::path& filePath) {
  using namespace MoleculeManip;

  std::cout << "Processing " << filePath.stem().string() << std::endl;

  // Read the file
  IO::MOLFileHandler molHandler;
  auto mol = molHandler.readSingle(filePath.string());

  CycleData cycleData {mol.getGraph()};

  testCycles(cycleData);

  auto indexMap = makeSmallestCycleMap(
    cycleData,
    mol
  );

  testCycles(cycleData);
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
