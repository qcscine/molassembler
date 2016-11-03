#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "IO.h"

using namespace MoleculeManip;

BOOST_AUTO_TEST_CASE(distancesMatrix) {
  // TODO this desperately needs more and better tests
  
  IO::MOLFileHandler molHandler;
  try {
    Molecule mol = molHandler.readSingle("mol_files/2,2-dimethybutane.mol");

    std::cout << "Distances Matrix: " << std::endl;
    std::cout << mol.getAdjacencyList().distancesMatrix() << std::endl << std::endl;

    auto distanceBounds = mol.getDistanceBoundsMatrix();

    std::cout << "Distance bounds matrix: " << std::endl;
    std::cout << distanceBounds << std::endl;
  } catch(const std::exception& e) {
    std::cout << e.what() << std::endl;
  }

}
