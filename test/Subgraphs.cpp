/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Subgraphs.h"
#include "molassembler/Molecule.h"
#include "molassembler/IO.h"
#include "molassembler/IO/SmilesParser.h"

#include <iostream>
#include "molassembler/Temple/Stringify.h"

using namespace Scine;
using namespace Molassembler;

BOOST_AUTO_TEST_CASE(SubgraphBasic) {
  const Molecule neopentane = IO::read("isomorphisms/neopentane.mol");
  const Molecule methyl = IO::experimental::parseSmilesSingleMolecule("[CH3]");

  const auto mappings = subgraphs::maximum(methyl, neopentane);
  BOOST_CHECK_MESSAGE(
    mappings.size() == 4,
    "Expected four mappings total, got " << mappings.size() << " instead"
  );
}
