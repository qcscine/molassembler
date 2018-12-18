/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"

#include "molassembler/Conformers.h"
#include "molassembler/Molecule.h"
#include "molassembler/IO.h"

#include "Utils/Typenames.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std::string_literals;
using namespace Scine;
using namespace molassembler;

BOOST_AUTO_TEST_CASE(FixedPositionsWork) {
  auto checkPositions = [](
    const Scine::Utils::PositionCollection& positions,
    const std::vector<std::pair<AtomIndex, Scine::Utils::Position>>& fixedPositions
  ) -> bool {
    bool pass = true;
    for(const auto& fixedPositionPair : fixedPositions) {
      if(
        !positions.row(fixedPositionPair.first).isApprox(
          fixedPositionPair.second,
          1e-4
        )
      ) {
        pass = false;
        std::cout << "Fixed position atom " << fixedPositionPair.first
          << " is at " << positions.row(fixedPositionPair.first).transpose()
          << ", but was supposed to be fixed at "
          << fixedPositionPair.second.transpose()
          << "\n";
      }
    }

    return pass;
  };

  auto octadecane = IO::read("various/octadecane.mol");

  // Start simple: make an arbitrary atom the origin.
  const Scine::Utils::Position origin(0.0, 0.0, 0.0);

  DistanceGeometry::Configuration config;
  config.fixedPositions = {{13, origin}};

  auto conformerResult = generateConformation(octadecane, config);
  BOOST_CHECK_MESSAGE(
    conformerResult,
    "Could not generate a conformer for octadecane with an atom fixed to the origin"
  );

  BOOST_CHECK_MESSAGE(
    checkPositions(conformerResult.value(), config.fixedPositions),
    "The fixed atom isn't approximately placed at the origin"
  );

  /* Octadecane carbon positions at both ends of the chain are 16, 17 (0-based).
   * Let's make this look ring-closing and force the rest of the dihedrals to
   * adapt.
   */
  config.fixedPositions = {
    {16, Scine::Utils::Position {-3, 0.0, 0.0}},
    {17, Scine::Utils::Position {3, 0.0, 0.0}}
  };
  conformerResult = generateConformation(octadecane, config);
  BOOST_CHECK_MESSAGE(
    conformerResult,
    "Could not generate a conformer for octadecane with ends close together"
  );

  BOOST_CHECK_MESSAGE(
    checkPositions(conformerResult.value(), config.fixedPositions),
    "The ring-like positions aren't fixed as required."
  );
}
