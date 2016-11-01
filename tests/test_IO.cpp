#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <iostream>
#include "IO.h"

using namespace MoleculeManip;

BOOST_AUTO_TEST_CASE( read_mol ) {
  // instantiate reader
  IO::MOLFileHandler molHandler;
  try {
    Molecule mol = molHandler.readSingle("mol_files/2,2-dimethybutane.mol");
    std::cout << mol << std::endl;
  } catch(const std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  try {
    Molecule mol = molHandler.readSingle("mol_files/asymCarbon.mol");
    std::cout << mol << std::endl;
  } catch(const std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
