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

BOOST_AUTO_TEST_CASE(SubgraphNeopentane, *boost::unit_test::label("Molassembler")) {
  const Molecule neopentane = IO::Experimental::parseSmilesSingleMolecule("CC(C)(C)C");
  const Molecule methyl = IO::Experimental::parseSmilesSingleMolecule("[CH3]");

  const auto maximumMappings = Subgraphs::maximum(methyl, neopentane);
  BOOST_CHECK_EQUAL(maximumMappings.size(), 4);
  for(const auto& mapping : maximumMappings) {
    // Each mapping should map four atoms (i.e. full size of the methyl)
    BOOST_CHECK_EQUAL(mapping.size(), 4);
  }

  const auto completeMappings = Subgraphs::complete(methyl, neopentane);
  BOOST_CHECK_EQUAL(completeMappings.size(), 4);
  for(const auto& mapping : completeMappings) {
    BOOST_CHECK_EQUAL(mapping.size(), 4);
  }
}

BOOST_AUTO_TEST_CASE(SubgraphOligopeptide, *boost::unit_test::label("Molassembler")) {
  const auto tetrapeptideSmiles = "CC(C)C(N)C(=O)NCC(=O)NC(CO)C(=O)NC(C)C(=O)O";
  const Molecule tetrapeptide = IO::Experimental::parseSmilesSingleMolecule(tetrapeptideSmiles);
  const auto bondPattern = IO::Experimental::parseSmilesSingleMolecule("[C](=O)[NH]");

  const auto completeMappings = Subgraphs::complete(bondPattern, tetrapeptide);
  // There should be three peptide bonds in a tetrapeptide
  BOOST_CHECK_EQUAL(completeMappings.size(), 3);
  for(const auto& mapping : completeMappings) {
    BOOST_CHECK_EQUAL(mapping.size(), 4);
  }

  // Arguments are not symmetric!
  BOOST_CHECK_EQUAL(Subgraphs::complete(tetrapeptide, bondPattern).size(), 0);
}
