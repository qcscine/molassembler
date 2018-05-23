#include "boost/program_options.hpp"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "DistanceGeometry/generateConformation.h"
#include "BoundsFromSymmetry.h"
#include "IO.h"
#include "AnalysisHelpers.h"
#include "StdlibTypeAlgorithms.h"
#include "Log.h"

#include "temple/constexpr/Numeric.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace molassembler;
using namespace molassembler::DistanceGeometry;

const std::string partialityChoices =
  "  0 - Four-Atom Metrization\n"
  "  1 - 10% Metrization\n"
  "  2 - All (default)\n";

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
    ("f", boost::program_options::value<std::string>(), "Read molecule to generate from file")
    ("i", boost::program_options::value<bool>(), "Specify whether inversion trick is to be used (Default: false)")
    ("p", boost::program_options::value<unsigned>(), "Set metrization partiality option (Default: full)")
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

    symmetries = {{Symmetry::allNames[argSymmetry]}};
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

  Partiality metrizationOption = Partiality::All;
  if(options_variables_map.count("p")) {
    unsigned index =  options_variables_map["p"].as<unsigned>();

    if(index > 2) {
      std::cout << "Specified metrization option is out of bounds. Valid choices are:\n"
        << partialityChoices;
      return 0;
    }

    metrizationOption = static_cast<DistanceGeometry::Partiality>(index);
  }

  Log::particulars.insert(Log::Particulars::DGStructureAcceptanceFailures);

/* Generating work */
  // Generate from file
  if(options_variables_map.count("f") == 1) {
    unsigned conformations = 1;
    if(options_variables_map.count("n")) {
      conformations = options_variables_map["n"].as<unsigned>();
    }

    bool useYInversionTrick = false;

    if(options_variables_map.count("i")) {
      useYInversionTrick = options_variables_map["i"].as<bool>();
    }

    auto filename = options_variables_map["f"].as<std::string>();

    if(!boost::filesystem::exists(filename)) {
      std::cout << "The specified file could not be found!" << std::endl;
      return 0;
    }

    auto mol = IO::read(filename);

    std::cout << mol << std::endl;

    auto debugData = detail::debugDistanceGeometry(
      mol,
      conformations,
      metrizationOption,
      useYInversionTrick
    );

    boost::filesystem::path filepath {filename};
    std::string filestem = filepath.stem().string();

    for(const auto& enumPair : enumerate(debugData)) {
      const auto& structNum = enumPair.index;
      const auto& refinementData = enumPair.value;

      std::string baseName = filestem + "-"s + std::to_string(structNum);

      AnalysisHelpers::writeDGPOVandProgressFiles(
        mol,
        baseName,
        refinementData
      );

      IO::write(
        filestem + "-"s + std::to_string(structNum) + "-last.mol"s,
        mol,
        DistanceGeometry::detail::convertToPositionCollection(
          refinementData.steps.back().positions
        )
      );
    }

    auto failures = temple::sum(
      temple::map(
        debugData,
        [](const auto& refinementData) -> unsigned {
          return static_cast<unsigned>(refinementData.isFailure);
        }
      )
    );

    if(failures > 0) {
      std::cout << "WARNING: " << failures << " refinements failed." << std::endl;
    }
  }

  // Not from file, then a basic molecule from symmetry
  if(options_variables_map.count("f") == 0) {

    for(const auto& symmetryName : symmetries) {

      // Make a molecule and generate an ensemble
      auto mol = DGDBM::asymmetricMolecule(symmetryName);

      auto debugData = detail::debugDistanceGeometry(
        mol,
        nStructures,
        metrizationOption,
        false,
        MoleculeSpatialModel::DistanceMethod::Uniform
      );

      for(const auto& enumPair : enumerate(debugData)) {
        const auto& structNum = enumPair.index;
        const auto& refinementData = enumPair.value;

        AnalysisHelpers::writeDGPOVandProgressFiles(
          mol,
          symmetryName,
          structNum,
          refinementData
        );

        IO::write(
          Symmetry::spaceFreeName(symmetryName) + "-"s
            + std::to_string(structNum) + "-last.mol"s,
          mol,
          DistanceGeometry::detail::convertToPositionCollection(
            refinementData.steps.back().positions
          )
        );
      }

      auto failures = temple::sum(
        temple::map(
          debugData,
          [](const auto& refinementData) -> unsigned {
            return static_cast<unsigned>(refinementData.isFailure);
          }
        )
      );

      if(failures > 0) {
        std::cout << "WARNING: " << failures << " refinements failed." << std::endl;
      }
    }
  }
}
