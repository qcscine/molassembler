#include "boost/program_options.hpp"

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"

#include "temple/constexpr/Numeric.h"
#include "temple/Containers.h"

#include "molassembler/Detail/AnalysisHelpers.h"
#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/DistanceGeometry/ConformerGeneration.h"
#include "molassembler/IO.h"
#include "molassembler/Log.h"

#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace molassembler;

const std::string partialityChoices =
  "  0 - Four-Atom Metrization\n"
  "  1 - 10% Metrization\n"
  "  2 - All (default)\n";

int main(int argc, char* argv[]) {
/* Set program options from command-line arguments */
  // Defaults
  unsigned nStructures = 1;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
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
  if(options_variables_map.count("help") > 0) {
    std::cout << options_description << std::endl;
    return 0;
  }

  if(options_variables_map.count("n") > 0) {
    unsigned argN = options_variables_map["n"].as<unsigned>();
    if(argN == 0) {
      std::cout << "Specified to generate zero structures. Exiting."
        << std::endl;
      return 0;
    }

    nStructures = argN;
  }

  DistanceGeometry::Partiality metrizationOption = DistanceGeometry::Partiality::All;
  if(options_variables_map.count("p") > 0) {
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
    bool useYInversionTrick = false;

    if(options_variables_map.count("i") > 0) {
      useYInversionTrick = options_variables_map["i"].as<bool>();
    }

    auto filename = options_variables_map["f"].as<std::string>();

    if(!boost::filesystem::exists(filename)) {
      std::cout << "The specified file could not be found!" << std::endl;
      return 0;
    }

    auto mol = IO::read(filename);

    std::cout << mol << std::endl;

    auto debugData = DistanceGeometry::debug(
      mol,
      nStructures,
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
        DistanceGeometry::detail::convertToAngstromWrapper(
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
