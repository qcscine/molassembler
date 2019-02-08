/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Editing.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"

BOOST_AUTO_TEST_CASE(cleaveEdit) {
  using namespace Scine;
  using namespace molassembler;

  if(!IO::LineNotation::enabled()) {
    BOOST_TEST_MESSAGE("obabel is not found. cleaveEdit is not run.");
    return;
  }

  Molecule caffeine = IO::LineNotation::fromCanonicalSMILES("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
  caffeine.canonicalize();

  const BondIndex bridge {15, 23};
  BOOST_REQUIRE_MESSAGE(!caffeine.graph().canRemove(bridge), "Selected edge for test is not a bridge edge");

  auto cleaved = Editing::cleave(caffeine, bridge);

  BOOST_CHECK_MESSAGE(
    cleaved.first.graph().N() + cleaved.second.graph().N() == caffeine.graph().N(),
    "Sum of number of vertices of cleaved molecules does not match that of the original molecule"
  );

  BOOST_CHECK_MESSAGE(
    cleaved.first.graph().B() + cleaved.second.graph().B() == caffeine.graph().B() - 1,
    "Sum of number of bonds of cleaved molecules does not match that of the original molecule minus one"
  );
}

BOOST_AUTO_TEST_CASE(insertEdit) {
  using namespace Scine;
  using namespace molassembler;

  if(!IO::LineNotation::enabled()) {
    BOOST_TEST_MESSAGE("obabel is not found. insertEdit is not run.");
    return;
  }

  // Set up the test and make sure the canonical indices still match expectations
  Molecule biphenyl = IO::LineNotation::fromCanonicalSMILES("C1=CC=C(C=C1)C2=CC=CC=C2");
  biphenyl.canonicalize();

  const BondIndex bridge {12, 13};
  BOOST_REQUIRE_MESSAGE(!biphenyl.graph().canRemove(bridge), "Selected edge for test is not a bridge edge");
  auto sides = biphenyl.graph().splitAlongBridge(bridge);

  BOOST_REQUIRE_MESSAGE(
    sides.first.size() == sides.second.size(),
    "Selected bridge edge does not split biphenyl equally in two"
  );

  Molecule pyrimidine = IO::LineNotation::fromCanonicalSMILES("C1=CN=CN=C1");
  pyrimidine.canonicalize();

  const AtomIndex firstPyrimidineNitrogen = 4, secondPyrimidineNitrogen = 5;
  BOOST_REQUIRE(pyrimidine.graph().elementType(firstPyrimidineNitrogen) == Utils::ElementType::N);
  BOOST_REQUIRE(pyrimidine.graph().elementType(secondPyrimidineNitrogen) == Utils::ElementType::N);

  // Insert pyrimidine into the bridge bond of biphenyl
  Molecule inserted = Editing::insert(
    biphenyl,
    pyrimidine,
    bridge,
    firstPyrimidineNitrogen,
    secondPyrimidineNitrogen
  );

  BOOST_CHECK_MESSAGE(
    inserted.graph().N() == biphenyl.graph().N() + pyrimidine.graph().N(),
    "Sum of number of vertices of log and wedge does not equal that of result"
  );

  BOOST_CHECK_MESSAGE(
    inserted.graph().B() == biphenyl.graph().B() + pyrimidine.graph().B() + 1,
    "Sum of number of edges of log and wedge plus one does not equal that of result"
  );
}

BOOST_AUTO_TEST_CASE(superposeEdit) {
  using namespace Scine;
  using namespace molassembler;

  if(!IO::LineNotation::enabled()) {
    BOOST_TEST_MESSAGE("obabel is not found. superposeEdit is not run.");
    return;
  }

  Molecule pyridine = IO::LineNotation::fromCanonicalSMILES("C1=CC=NC=C1");
  pyridine.canonicalize();
  const AtomIndex pyridineNitrogen = 5;
  BOOST_REQUIRE(pyridine.graph().elementType(pyridineNitrogen) == Utils::ElementType::N);

  Molecule methane = IO::LineNotation::fromCanonicalSMILES("C");
  methane.canonicalize();
  const AtomIndex methaneHydrogen = 0;
  BOOST_REQUIRE(methane.graph().elementType(methaneHydrogen) == Utils::ElementType::H);

  Molecule superposition = Editing::superpose(
    pyridine,
    methane,
    pyridineNitrogen,
    methaneHydrogen
  );

  BOOST_CHECK(superposition.graph().N() == pyridine.graph().N() + methane.graph().N() - 1);
}

BOOST_AUTO_TEST_CASE(substituteEdit) {
  using namespace Scine;
  using namespace molassembler;

  if(!IO::LineNotation::enabled()) {
    BOOST_TEST_MESSAGE("obabel is not found. substituteEdit is not run.");
    return;
  }

  Molecule chlorobenzene = IO::LineNotation::fromCanonicalSMILES("C1=CC=C(C=C1)Cl");
  chlorobenzene.canonicalize();
  const BondIndex chlorobenzeneBridge {5, 7};
  BOOST_REQUIRE_MESSAGE(!chlorobenzene.graph().canRemove(chlorobenzeneBridge), "Selected edge for test is not a bridge edge");

  Molecule phenole = IO::LineNotation::fromCanonicalSMILES("C1=CC=C(C=C1)O");
  phenole.canonicalize();
  const BondIndex phenoleBridge {6, 8};
  BOOST_REQUIRE_MESSAGE(!phenole.graph().canRemove(phenoleBridge), "Selected edge for test is not a bridge edge");

  Molecule substituted = Editing::substitute(
    chlorobenzene,
    phenole,
    chlorobenzeneBridge,
    phenoleBridge
  );

  Molecule biphenyl = IO::LineNotation::fromCanonicalSMILES("C1=CC=C(C=C1)C2=CC=CC=C2");
  BOOST_CHECK(substituted == biphenyl);

  substituted.canonicalize();
  biphenyl.canonicalize();
  BOOST_CHECK(substituted == biphenyl);
}

BOOST_AUTO_TEST_CASE(connectEdit) {
  using namespace Scine;
  using namespace molassembler;

  if(!IO::LineNotation::enabled()) {
    BOOST_TEST_MESSAGE("obabel is not found. connectEdit is not run.");
    return;
  }

  Molecule pyridine = IO::LineNotation::fromCanonicalSMILES("C1=CC=NC=C1");
  pyridine.canonicalize();
  const AtomIndex pyridineNitrogen = 5;
  BOOST_REQUIRE(pyridine.graph().elementType(pyridineNitrogen) == Utils::ElementType::N);

  Molecule connected = Editing::connect(
    pyridine,
    pyridine,
    pyridineNitrogen,
    pyridineNitrogen,
    BondType::Single
  );

  BOOST_CHECK(connected.graph().N() == 2 * pyridine.graph().N());
}
