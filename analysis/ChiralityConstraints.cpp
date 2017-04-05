#include "IO.h"
#include "Log.h"
#include <iostream>

/* This analysis requires a full set of DGRefinement analysis data. It will
 * overwrite DGRefinement's <symmetry>.csv files with the symmetry fit log data
 * from SymmetryFit's operator << for all generated optimized files
 */

int main() {
  // Set the define to log full stereocenter fitting information, and nothing else
  Log::level = Log::Level::None;
  Log::particulars = {Log::Particulars::StereocenterFitAnalysisInfo};

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
        + ".mol"s;

      Molecule mol = molHandler.readSingle(filename);
    }

    // Restore cout
    std::cout.rdbuf(coutBuf);
  }

  return 0;
}
