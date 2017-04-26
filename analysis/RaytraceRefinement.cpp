#include "boost/program_options.hpp"

#include "DistanceGeometry/generateConformation.h"
#include "BoundsFromSymmetry.h"
#include "IO.h"
#include "Log.h"
#include "template_magic/Enumerate.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace MoleculeManip;
using namespace MoleculeManip::DistanceGeometry;

inline char mapIndexToChar(const unsigned& index) {
  return 'A' + index;
}

int main(int argc, char* argv[]) {
/* Options not altered by command-line arguments */
  const unsigned dimensionality = 4;

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
    // Make a space-free string from the name
    std::string spaceFreeName = Symmetry::name(symmetryName);
    std::replace(
      spaceFreeName.begin(),
      spaceFreeName.end(),
      ' ',
      '-'
    );

    // Make a molecule and generate an ensemble
    auto mol = DGDBM::symmetricMolecule(symmetryName);

    auto debugData = detail::debugDistanceGeometry(
      mol,
      nStructures,
      MetrizationOption::off,
      false
    );

    // Write the POV-Ray files
    for(const auto& enumPair : enumerate(debugData.refinements)) {
      const auto& structNum = enumPair.index;
      const auto& refinementData = enumPair.value;

      std::string progressFilename = spaceFreeName + "-"s 
        + std::to_string(structNum) + "-progress.csv"s;
      std::ofstream progressFile (progressFilename);

      progressFile << std::scientific;

      for(const auto& refinementEnumPair : enumerate(refinementData.steps)) {

        const auto& stepData = refinementEnumPair.value;
        assert(stepData.positions.size() % dimensionality == 0);
        const unsigned N = stepData.positions.size() / dimensionality;

        // Collect sum of absolute 4D values
        double totalAbs4D = 0;
        for(unsigned i = 0; i < N; i++) {
          totalAbs4D += std::fabs(stepData.positions[4 * i + 3]);
        }

        progressFile << stepData.error << "," 
          << stepData.gradient.norm() << "," 
          << static_cast<unsigned>(stepData.compress) << "," 
          << totalAbs4D << "\n";

        std::stringstream filename;
        filename << spaceFreeName << "-"
          << structNum << "-"
          << std::setfill('0') << std::setw(3) << refinementEnumPair.index
          << ".pov";

        std::ofstream outStream(
          filename.str()
        );

        outStream << "#version 3.7;\n"
          << "#include \"scene.inc\"\n\n";

        // Define atom names with positions
        for(unsigned i = 0; i < N; i++) {
          outStream << "#declare " << mapIndexToChar(i) <<  " = <";
          outStream << std::fixed << std::setprecision(4);
          outStream << stepData.positions[dimensionality * i] << ", ";
          outStream << stepData.positions[dimensionality * i + 1] << ", ";
          outStream << stepData.positions[dimensionality * i + 2] << ">;\n";
        }
        outStream << "\n";

        // Atoms
        for(unsigned i = 0; i < N; i++) {
          outStream << "Atom4D(" << mapIndexToChar(i) << ", " 
            << std::fabs(stepData.positions[dimensionality * i + 3]) << ")\n";
        }
        outStream << "\n";

        // Bonds
        for(const auto& edgePair : mol.getAdjacencyList().getEdges()) {
          outStream << "Bond(" 
            << mapIndexToChar(edgePair.first.first) << ","
            << mapIndexToChar(edgePair.first.second) 
            << ")\n";
        }
        outStream << "\n";

        // Tetrahedra
        if(refinementData.constraints.size() > 0) {
          for(const auto& chiralityConstraint : refinementData.constraints) {
            outStream << "TetrahedronHighlight("
              << mapIndexToChar(chiralityConstraint.indices[0]) << ", "
              << mapIndexToChar(chiralityConstraint.indices[1]) << ", "
              << mapIndexToChar(chiralityConstraint.indices[2]) << ", "
              << mapIndexToChar(chiralityConstraint.indices[3]) 
              << ")\n";
          }
          outStream << "\n";
        }

        // Gradients
        for(unsigned i = 0; i < N; i++) {
          outStream << "GradientVector(" << mapIndexToChar(i) <<", <"
            << std::fixed << std::setprecision(4)
            << (-stepData.gradient[3 * i]) << ", "
            << (-stepData.gradient[3 * i + 1]) << ", "
            << (-stepData.gradient[3 * i + 2]) << ", "
            << ">)\n";
        }

        outStream.close();

      }
    }
  }
}
