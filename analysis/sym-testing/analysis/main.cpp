#include "MoleculeAlgorithms.h"
#include "IO.h"

#include <iostream>

int main() {
  using namespace MoleculeManip;
  using namespace std::string_literals;

  // instantiate reader
  IO::MOLFileHandler molHandler;

  for(const auto& symmetryName : Symmetry::allNames) {
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    // Open a csv file for the symmetry
    std::ofstream outFile(spaceFreeName+".csv"s);
    // Save so we can restore
    std::streambuf* coutBuf = std::cout.rdbuf();

    // Redirect cout to outFile
    std::cout.rdbuf(outFile.rdbuf());


    for(unsigned i = 0; i < 100; i++) {
      std::string filename = "opt-"s +spaceFreeName + std::to_string(i) 
        + ".xyz"s;

      try {
        Molecule mol = molHandler.readSingle(filename);
      } catch(const std::exception& e) {
        std::cout << e.what() << std::endl;
        throw e;
      }

    }

    // Restore cout
    std::cout.rdbuf(coutBuf);
  }

  return 0;
}
