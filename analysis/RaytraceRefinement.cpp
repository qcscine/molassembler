#include "boost/program_options.hpp"

#include "DistanceGeometry/generateConformation.h"
#include "BoundsFromSymmetry.h"
#include "IO.h"
#include "Log.h"
#include "AnalysisHelpers.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

int main(int argc, char* argv[]) {
/* Set program options from command-line arguments */
  // Defaults
  unsigned nStructures = 1;
  auto symmetries = Symmetry::allNames;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    ("s", boost::program_options::value<unsigned>(), "Specify symmetry index (zero-based)")
    ("n", boost::program_options::value<unsigned>(), "Set number of structures to generate")
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, options_description),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);  

  // Manage the results
  if(options_variables_map.count("help")) {
    std::cout << options_description << std::endl;
    return 0;
  }

  if(options_variables_map.count("s")) {
    unsigned argSymmetry = options_variables_map["s"].as<unsigned>();
    if(argSymmetry >= Symmetry::allNames.size()) {
      std::cout << "Specified symmetry out of bounds. Valid symmetries are 0-"
        << (Symmetry::allNames.size() - 1) << ":\n\n";
      for(unsigned i = 0; i < Symmetry::allNames.size(); i++) {
        std::cout << "  " << i << " - " << Symmetry::name(Symmetry::allNames.at(i)) << "\n";
      }
      std::cout << std::endl;
      return 0;
    }

    symmetries = {Symmetry::allNames[argSymmetry]};
  }

  if(options_variables_map.count("n")) {
    unsigned argN = options_variables_map["n"].as<unsigned>();
    if(argN == 0) {
      std::cout << "Specified to generate zero structures. Exiting."
        << std::endl;
      return 0;
    }

    nStructures = argN;
  }

/* Generating work */
  Log::level = Log::Level::None;
  Log::particulars = {Log::Particulars::DGRefinementChiralityNumericalDebugInfo};

  for(const auto& symmetryName : symmetries) {

    // Make a molecule and generate an ensemble
    auto mol = DGDBM::symmetricMolecule(symmetryName);

    auto debugData = detail::debugDistanceGeometry(
      mol,
      nStructures,
      MetrizationOption::full,
      false
    );

    for(const auto& enumPair : enumerate(debugData.refinements)) {
      const auto& structNum = enumPair.index;
      const auto& refinementData = enumPair.value;

      AnalysisHelpers::writeDGPOVandProgressFiles(
        mol,
        symmetryName,
        structNum,
        refinementData
      );
    }

  }
}
