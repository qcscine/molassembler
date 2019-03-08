/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "molassembler/DirectedConformerGenerator.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule.h"
#include "temple/Stringify.h"

using namespace std::string_literals;
using namespace Scine;
using namespace molassembler;

const std::string partialityChoices =
  "  0 - Four-Atom Metrization (default)\n"
  "  1 - 10% Metrization\n"
  "  2 - All\n";

int main(int argc, char* argv[]) {
  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help,h", "Produce help message")
    (
      "file,f",
      boost::program_options::value<std::string>(),
      "Read molecule to generate from file"
    )
    (
      "partiality,p",
      boost::program_options::value<unsigned>(),
      "Set metrization partiality option (Default: four-atom)"
    )
    (
      "steps,s",
      boost::program_options::value<unsigned>(),
      "Alter the maximum number of refinement steps (Default: 10'000)"
    )
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::positional_options_description positional_description;
  positional_description.add("file", 1);
  boost::program_options::store(
    boost::program_options::command_line_parser(argc, argv).
    options(options_description).
    positional(positional_description).
    style(
      boost::program_options::command_line_style::unix_style
      | boost::program_options::command_line_style::allow_long_disguise
    ).run(),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  if(options_variables_map.count("help") > 0 || options_variables_map.count("file") == 0) {
    std::cout << options_description << std::endl;
    return 0;
  }

  DistanceGeometry::Configuration configuration;
  if(options_variables_map.count("partiality") > 0) {
    unsigned index =  options_variables_map["partiality"].as<unsigned>();

    if(index > 2) {
      std::cout << "Specified metrization option is out of bounds. Valid choices are:\n"
        << partialityChoices;
      return 0;
    }

    configuration.partiality = static_cast<DistanceGeometry::Partiality>(index);
  }

  if(options_variables_map.count("steps") > 0) {
    configuration.refinementStepLimit = options_variables_map["steps"].as<unsigned>();
  }

  std::string filename = options_variables_map["file"].as<std::string>();

  Molecule mol = IO::read(filename);

  std::cout << mol << "\n";

  boost::filesystem::path filepath {filename};
  std::string filestem = filepath.stem().string();

  DirectedConformerGenerator generator(mol);

  if(generator.bondList().empty()) {
    std::cout << "No bonds in this molecule are relevant to directed conformer generation.\n";

    const unsigned maxTries = 3;
    for(unsigned attempt = 0; attempt < maxTries; ++attempt) {
      auto positionResult = generateConformation(mol);
      if(positionResult) {
        std::cout << "Generated conformation.\n";
        IO::write(
          filestem + "-0.mol",
          mol,
          positionResult.value()
        );
        break;
      } else {
        std::cout << "Could not generate conformer: "
          << positionResult.error().message() << "\n";
      }
    }

    return 0;
  }

  for(const auto& bondIndex : generator.bondList()) {
    std::cout << "Considering bond " << bondIndex.first << "-" << bondIndex.second << "\n";
  }

  std::cout << "Ideally need " << generator.idealEnsembleSize() << " conformers.\n";

  const unsigned maxTries = 3;

  unsigned conformerCount = 0;
  while(generator.decisionListSetSize() != generator.idealEnsembleSize()) {
    auto newDecisionList = generator.generateNewDecisionList();

    for(unsigned attempt = 0; attempt < maxTries; ++attempt) {
      const auto& confMol = generator.conformationMolecule(newDecisionList);
      auto positionResult = generateConformation(confMol, configuration);

      if(positionResult) {
        std::cout << "Generated conformer #" << (conformerCount + 1) << ", decision list " << temple::stringify(newDecisionList) << "\n";
        IO::write(
          filestem + "-"s + std::to_string(conformerCount + 1) + ".mol",
          mol,
          positionResult.value()
        );
        break;
      }

      std::cout << "Could not generate decision list " << temple::stringify(newDecisionList) << ": " << positionResult.error().message() << "\n";
      IO::write(
        filestem + "-"s + std::to_string(conformerCount + 1) + ".masm",
        confMol
      );
    }

    ++conformerCount;
  }

  return 0;
}
