/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/test/unit_test.hpp"

#include "Molassembler/Conformers.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/IO.h"

#include "Utils/Typenames.h"

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std::string_literals;
using namespace Scine;
using namespace Molassembler;

BOOST_AUTO_TEST_CASE(FixedPositionsWork, *boost::unit_test::label("DG")) {
  auto checkPositions = [](
    const Utils::PositionCollection& positions,
    const std::vector<std::pair<AtomIndex, Utils::Position>>& fixedPositions
  ) -> bool {
    bool pass = true;
    for(const auto& fixedPositionPair : fixedPositions) {
      if(
        !positions.row(fixedPositionPair.first).isApprox(
          fixedPositionPair.second,
          1e-2
        )
      ) {
        pass = false;
        std::cout << "Fixed position atom " << fixedPositionPair.first
          << " is at " << positions.row(fixedPositionPair.first)
          << ", but was supposed to be fixed at "
          << fixedPositionPair.second
          << "\n";
      }
    }

    return pass;
  };

  auto octadecane = IO::read("various/octadecane.mol");

  // Start simple: make an arbitrary atom the origin.
  const Utils::Position origin(0.0, 0.0, 0.0);

  DistanceGeometry::Configuration config;
  config.fixedPositions = {{13, origin}};

  auto conformerResult = generateRandomConformation(octadecane, config);
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
    {16, Utils::Position {-3, 0.0, 0.0}},
    {17, Utils::Position {3, 0.0, 0.0}}
  };
  conformerResult = generateRandomConformation(octadecane, config);
  BOOST_CHECK_MESSAGE(
    conformerResult,
    "Could not generate a conformer for octadecane with ends close together"
  );

  BOOST_CHECK_MESSAGE(
    checkPositions(conformerResult.value(), config.fixedPositions),
    "The ring-like positions aren't fixed as required."
  );
}
