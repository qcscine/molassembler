#define BOOST_TEST_MODULE MoleculeTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include "IO.h"

BOOST_AUTO_TEST_CASE( read_mol ) {
  using namespace MoleculeManip;

  std::vector<std::string> files {
    "../tests/mol_files/2,2-dimethybutane.mol",
    "../tests/mol_files/asymCarbon.mol",
    "../tests/mol_files/C8H12_asym.mol",
  };

  // instantiate reader
  IO::MOLFileHandler molHandler;

  for(const auto& filename : files) {
    try {
      Molecule mol = molHandler.readSingle(filename);
      // Invoke ostream operator
      std::cout << mol << std::endl;

      // Make dot files for every file
      auto slashSplat = StdlibTypeAlgorithms::split(filename, '/');
      auto dotSplat = StdlibTypeAlgorithms::split(slashSplat.back(), '.');

      mol.dumpGraphviz(dotSplat.front() + ".dot");

    } catch(const std::exception& e) {
      std::cout << e.what() << std::endl;
      throw e;
    }
  }
}
