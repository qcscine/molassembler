/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <algorithm>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <functional>

#include "boost/program_options.hpp"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Shapes/Data.h"

#include "Molassembler/Stereopermutation/Manipulation.h"

std::ostream& nl(std::ostream& os) {
  os << '\n';
  return os;
}

using namespace Scine;
using namespace Molassembler;
using namespace Stereopermutations;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
  // Set up option parsing
  boost::program_options::options_description options_description("Recognized options");
  options_description.add_options()
    ("help", "Produce help message")
    (
      "shape,s",
      boost::program_options::value<unsigned>(),
      "Set shape index"
    )
    (
      "characters,c",
      boost::program_options::value<std::string>(),
      "Specify case characters"
    )
    (
      "links,l",
      boost::program_options::value<std::string>(),
      "Specify links (remember to place quotation marks)"
    )
  ;

  // Parse
  boost::program_options::variables_map options_variables_map;
  boost::program_options::positional_options_description positional_description;
  positional_description.add("shape", 1);
  positional_description.add("characters", 1);
  positional_description.add("links", 1);

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

  if(options_variables_map.count("help") != 0) {
    std::cout << options_description << nl;
    std::cout << "Example flags for octahedral (A-A)3: -s 12 -c AAAAAA -l \"0, 1, 2, 3, 4, 5\"" << nl;
    std::cout << "Note: The links list is stripped of anything not a comma or number,  split by commas and parsed pairwise." << nl;

    std::cout << nl << "Valid shape indices:" << nl << nl;
    for(unsigned i = 0; i < Shapes::allShapes.size(); i++) {
      std::cout << std::setw(4) << i << " - " << Shapes::name(Shapes::allShapes.at(i)) << nl;
    }
    std::cout << nl;
    return 0;
  }

  if(options_variables_map.count("shape") > 0 && options_variables_map.count("characters") > 0) {
    // Validate shape argument
    unsigned argSymmetry = options_variables_map["shape"].as<unsigned>();
    if(argSymmetry >= Shapes::allShapes.size()) {
      std::cout << "Specified shape out of bounds. Valid shapes are 0-"
        << (Shapes::allShapes.size() - 1) << ":" << nl << nl;
      for(unsigned i = 0; i < Shapes::allShapes.size(); i++) {
        std::cout << "  " << i << " - " << Shapes::name(Shapes::allShapes.at(i)) << nl;
      }
      std::cout << nl;
      return 0;
    }

    Shapes::Shape shape = Shapes::allShapes[argSymmetry];

    // Validate characters
    std::string chars = options_variables_map["characters"].as<std::string>();

    if(chars.size() != Shapes::size(shape)) {
      std::cout << "Number of characters does not fit specified shape size: "
        << chars.size() << " characters specified, shape size is "
        << Shapes::size(shape) << nl;

      return 0;
    }

    // Validate links (if present)
    Stereopermutation::OrderedLinks links;
    if(options_variables_map.count("links") != 0) {
      /* Naive parse strategy:
       * - remove anything not a number or comma from the string
       * - split along commas, check size % 2 == 0
       * - parse non-overlapping pairwise as links
       */
      std::string linksString = options_variables_map["links"].as<std::string>();

      linksString.erase(
        std::remove_if(
          linksString.begin(),
          linksString.end(),
          [](const char& a) -> bool {
            return !(
              a == ',' || (std::isdigit(a) != 0)
            );
          }
        ),
        linksString.end()
      );

      std::vector<unsigned> linkIndices;

      std::stringstream ss;
      ss.str(linksString);
      std::string item;
      while(std::getline(ss, item, ',')) {
        linkIndices.push_back(std::stoul(item));
      }

      if(linkIndices.size() % 2 != 0) {
        std::cout << "The number of provided link indices is not even. You must "
          << "specify link pairs as pairwise indices for each link" << nl;
        return 0;
      }

      for(unsigned i = 0; i < linkIndices.size() / 2; ++i) {
        const auto& a = linkIndices.at(2 * i);
        const auto& b = linkIndices.at(2 * i + 1);

        if(a == b) {
          std::cout << "You have specified a link between identical indices!" << nl;
          return 0;
        }

        links.emplace_back(
          std::min(a, b),
          std::max(a, b)
        );
      }
    }

    // Generate the assignment
    Stereopermutation base {Stereopermutation::occupationFromChars(chars), links};

    auto unique = uniques(base, shape, false);

    std::cout << "Shape: " << Shapes::name(shape) << nl
      << "Characters: " << chars << nl;

    if(!links.empty()) {
      std::cout << "Links: " << Temple::stringifyContainer(
        links,
        [](const Stereopermutation::Link link) -> std::string {
          return "{"s + std::to_string(link.first) + ", "s + std::to_string(link.second) + "}"s;
        }
      ) << nl;
    }

    std::cout << nl;

    for(unsigned i = 0; i < unique.list.size(); ++i) {
      std::cout << "Weight " << unique.weights[i] << ": "
        << unique.list[i].toString();

      if(!unique.list[i].links.empty()) {
        std::cout << ", link angles: ";
      }

      for(const auto& linkPair : unique.list[i].links) {
        std::cout << (
          180 * Shapes::angleFunction(shape)(linkPair.first, linkPair.second) / M_PI
        ) << " ";
      }
      std::cout << nl;
    }

    std::cout << unique.list.size() << " stereopermutations\n";
  }

  return 0;
}
