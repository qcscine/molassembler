/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/IO/SmilesParser.h"
#include "Molassembler/IO.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Descriptors.h"

// #include "Molassembler/Temple/Stringify.h"

BOOST_AUTO_TEST_CASE(RankingEquivalenceTests, *boost::unit_test::label("Molassembler")) {
  using namespace Scine;
  using namespace Molassembler;

  const std::vector<std::pair<std::string, unsigned>> expectations {
    {"C", 2},
    {"CC", 4}, // NOTE: 2, but no bond symmetry
    {"CCC", 4},
    {"CCCC", 8}, // NOTE: 4, but no bond symmetry
    {"CCCCC", 6},
    {"C1CC1", 2},
    {"C1CCC1", 4}, // NOTE: 2, but no bond symmetry
    {"C1CCCC1", 2},
    {"C1CCCCC1", 4}, // NOTE: 2, but bond symmetry
    {"C[Fe]1[H][H]1", 4},  // H3C-Fe-(eta2-H2)
  };

  // unsigned count = 0;
  for(const auto& iterPair : expectations) {
    const auto molecule = IO::Experimental::parseSmilesSingleMolecule(iterPair.first);
    const auto nonEquivalent = rankingDistinctAtoms(molecule);
    BOOST_TEST_CONTEXT("For smiles: " << iterPair.first) {
      BOOST_CHECK_EQUAL(nonEquivalent.size(), iterPair.second);
    }

    // if(nonEquivalent.size() != iterPair.second) {
    //   IO::write("failure-" + std::to_string(count) + ".svg", molecule);
    //   std::cout << Temple::stringify(nonEquivalent) << "\n";
    //   ++count;
    // }
  }
}
