/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/GraphAlgorithms.h"
#include "Molassembler/IO/SmilesParser.h"
#include "Molassembler/IO.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Graph/EditDistance.h"

#include <iostream>
#include "Molassembler/Temple/Stringify.h"

using namespace Scine;
using namespace Molassembler;

BOOST_AUTO_TEST_CASE(PathsGeneration, *boost::unit_test::label("Molassembler")) {
  const std::string smiles = "C1CCCC1";
  const auto cyclopentane = IO::Experimental::parseSmilesSingleMolecule(smiles);

  BOOST_REQUIRE(
    Temple::all_of(
      std::vector<AtomIndex> {{0, 1, 2, 3, 4}},
      [&](const AtomIndex i) -> bool {
        return cyclopentane.graph().elementType(i) == Utils::ElementType::C;
      }
    )
  );

  const auto predecessors = shortestPaths(0, cyclopentane.graph());

  const auto toAdjacent = predecessors.path(1);
  BOOST_CHECK_EQUAL(toAdjacent.size(), 2);
  BOOST_CHECK_EQUAL(toAdjacent.front(), 0);
  BOOST_CHECK_EQUAL(toAdjacent.back(), 1);

  const auto toExtended = predecessors.path(2);
  BOOST_CHECK_EQUAL(toExtended.size(), 3);
  BOOST_CHECK_EQUAL(toExtended.at(0), 0);
  BOOST_CHECK_EQUAL(toExtended.at(1), 1);
  BOOST_CHECK_EQUAL(toExtended.at(2), 2);
}

BOOST_AUTO_TEST_CASE(UniqueDescendants, *boost::unit_test::label("Molassembler")) {
  const auto parse = &IO::Experimental::parseSmilesSingleMolecule;
  const auto cyclopentane = parse("C1CCCC1");
  const auto cyclohexane = parse("C1CCCCC1");

  const auto getCarbonNeighbors = [](const AtomIndex i, const Molecule& m) {
    std::vector<AtomIndex> vs;
    for(AtomIndex adjacent : m.graph().adjacents(i)) {
      if(m.graph().elementType(adjacent) == Utils::ElementType::C) {
        vs.push_back(adjacent);
      }
    }
    return vs;
  };

  const auto cyclopentaneNeighbors = getCarbonNeighbors(0, cyclopentane);
  const auto pentaneDescendants = GraphAlgorithms::bfsUniqueDescendants(
    0,
    cyclopentaneNeighbors,
    cyclopentane.graph().inner()
  );

  const auto count = [](const AtomIndex v, const std::vector<AtomIndex>& vs) {
    return std::count(std::begin(vs), std::end(vs), v);
  };

  BOOST_CHECK_EQUAL(count(0, pentaneDescendants), 3);
  BOOST_CHECK_EQUAL(count(cyclopentaneNeighbors.front(), pentaneDescendants), 6);
  BOOST_CHECK_EQUAL(count(cyclopentaneNeighbors.back(), pentaneDescendants), 6);

  const auto cyclohexaneNeighbors = getCarbonNeighbors(0, cyclohexane);
  const auto hexaneDescendants = GraphAlgorithms::bfsUniqueDescendants(
    0,
    cyclohexaneNeighbors,
    cyclohexane.graph().inner()
  );

  constexpr auto split = std::numeric_limits<AtomIndex>::max();
  BOOST_CHECK_EQUAL(count(0, hexaneDescendants), 3);
  BOOST_CHECK_EQUAL(count(split, hexaneDescendants), 3);
  BOOST_CHECK_EQUAL(count(cyclohexaneNeighbors.front(), hexaneDescendants), 6);
  BOOST_CHECK_EQUAL(count(cyclohexaneNeighbors.back(), hexaneDescendants), 6);
}

BOOST_AUTO_TEST_CASE(EditDistance, *boost::unit_test::label("Molassembler")) {
  const auto parse = &IO::Experimental::parseSmilesSingleMolecule;
  auto methane = parse("C");
  auto methyl = parse("[CH3]");
  BOOST_CHECK_EQUAL(
    editDistance(methane.graph(), methyl.graph()),
    2
  );
  // Must be symmetric
  BOOST_CHECK_EQUAL(
    editDistance(methyl.graph(), methane.graph()),
    2
  );

  auto silane = parse("[SiH4]");
  BOOST_CHECK_EQUAL(
    editDistance(silane.graph(), methane.graph()),
    1
  );
  BOOST_CHECK_EQUAL(
    editDistance(silane.graph(), methyl.graph()),
    3
  );
}
