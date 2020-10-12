/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Subgraphs.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/IO.h"
#include "Molassembler/IO/SmilesParser.h"

#include <iostream>
#include "Molassembler/Temple/Stringify.h"

using namespace Scine;
using namespace Molassembler;

BOOST_AUTO_TEST_CASE(SubgraphBasic, *boost::unit_test::label("Molassembler")) {
  const Molecule neopentane = IO::read("isomorphisms/neopentane.mol");
  const Molecule methyl = IO::Experimental::parseSmilesSingleMolecule("[CH3]");

  const auto mappings = subgraphs::maximum(methyl, neopentane);
  BOOST_CHECK_MESSAGE(
    mappings.size() == 4,
    "Expected four mappings total, got " << mappings.size() << " instead"
  );
}
