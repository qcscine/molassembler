/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"

#include "Molassembler/Molecule.h"
#include "Molassembler/IO.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"

#include <iostream>

int main(int argc, char* argv[]) {
  using namespace Scine;
  using namespace Molassembler;

  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    ("f", boost::program_options::value<std::string>(), "Which file to make an isomorphism to.")
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::store(
    boost::program_options::parse_command_line(argc, argv, options_description),
    options_variables_map
  );
  boost::program_options::notify(options_variables_map);

  if(options_variables_map.count("f") > 0) {
    boost::filesystem::path filepath {
      options_variables_map["f"].as<std::string>()
    };

    if(!boost::filesystem::exists(filepath)) {
      std::cout << "That file does not exist!" << std::endl;
      return 0;
    }

    // This can throw in lots of cases
    auto readData = Utils::ChemicalFileHandler::read(filepath.string());

    auto shuffledData = IO::shuffle(readData.first, readData.second);

    Utils::ChemicalFileHandler::write(
      filepath.stem().string() + "isomorphism_.mol",
      std::get<0>(shuffledData),
      std::get<1>(shuffledData)
    );
  } else {
    std::cout << options_description << std::endl;
  }

  return 0;
}
