#include "BoostTestingHeader.h"
#include <iostream>
#include "IO.h"

BOOST_AUTO_TEST_CASE( read_mol ) {
  using namespace MoleculeManip;

  // instantiate reader
  IO::MOLFileHandler molHandler;
  try {
    Molecule mol = molHandler.readSingle("../tests/mol_files/2,2-dimethybutane.mol");
    std::cout << mol << std::endl;
  } catch(const std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  try {
    Molecule mol = molHandler.readSingle("../tests/mol_files/asymCarbon.mol");
    std::cout << mol << std::endl;
  } catch(const std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  try {
    Molecule mol = molHandler.readSingle("../tests/mol_files/C8H12_asym.mol");
    std::cout << mol << std::endl;
  } catch(const std::exception& e) {
    std::cout << e.what() << std::endl;
  }
}
