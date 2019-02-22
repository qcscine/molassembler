/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>

#include "molassembler/Editing.h"
#include "molassembler/IO.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"

/* SMILES for molecules imported here
 * - caffeine: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
 * - biphenyl: "C1=CC=C(C=C1)C2=CC=CC=C2"
 * - pyrimidine: "C1=CN=CN=C1"
 * - pyridine:"C1=CC=NC=C1"
 * - methane: "C"
 * - chlorobenzene: "C1=CC=C(C=C1)Cl"
 * - phenole: "C1=CC=C(C=C1)O"
 * - mesitylene: "CC1=CC(=CC(=C1)C)C"
 * - nhc: "C1NC=CN1"
 *
 * The .masm are just IO::LineNotation::fromCanonicalSMILES and directly
 * exported using IO::write(str, mol).
 */

BOOST_AUTO_TEST_CASE(cleaveEdit) {
  using namespace Scine;
  using namespace molassembler;

  Molecule caffeine = IO::read("masm/caffeine.masm");
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

  // Set up the test and make sure the canonical indices still match expectations
  Molecule biphenyl = IO::read("masm/biphenyl.masm");
  biphenyl.canonicalize();

  const BondIndex bridge {12, 13};
  BOOST_REQUIRE_MESSAGE(!biphenyl.graph().canRemove(bridge), "Selected edge for test is not a bridge edge");
  auto sides = biphenyl.graph().splitAlongBridge(bridge);

  BOOST_REQUIRE_MESSAGE(
    sides.first.size() == sides.second.size(),
    "Selected bridge edge does not split biphenyl equally in two"
  );

  Molecule pyrimidine = IO::read("masm/pyrimidine.masm");
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

  Molecule pyridine = IO::read("masm/pyridine.masm");
  pyridine.canonicalize();
  const AtomIndex pyridineNitrogen = 5;
  BOOST_REQUIRE(pyridine.graph().elementType(pyridineNitrogen) == Utils::ElementType::N);

  Molecule methane = IO::read("masm/methane.masm");
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

  Molecule chlorobenzene = IO::read("masm/chlorobenzene.masm");
  chlorobenzene.canonicalize();
  const BondIndex chlorobenzeneBridge {5, 7};
  BOOST_REQUIRE_MESSAGE(!chlorobenzene.graph().canRemove(chlorobenzeneBridge), "Selected edge for test is not a bridge edge");

  Molecule phenole = IO::read("masm/phenole.masm");
  phenole.canonicalize();
  const BondIndex phenoleBridge {6, 8};
  BOOST_REQUIRE_MESSAGE(!phenole.graph().canRemove(phenoleBridge), "Selected edge for test is not a bridge edge");

  Molecule substituted = Editing::substitute(
    chlorobenzene,
    phenole,
    chlorobenzeneBridge,
    phenoleBridge
  );

  Molecule biphenyl = IO::read("masm/biphenyl.masm");
  BOOST_CHECK(substituted == biphenyl);

  substituted.canonicalize();
  biphenyl.canonicalize();
  BOOST_CHECK(substituted == biphenyl);
}

BOOST_AUTO_TEST_CASE(connectEdit) {
  using namespace Scine;
  using namespace molassembler;

  Molecule pyridine = IO::read("masm/pyridine.masm");
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

BOOST_AUTO_TEST_CASE(mesitylSubstitutionBug) {
  using namespace Scine;
  using namespace molassembler;

  Molecule mesitylene = IO::read("masm/mesitylene.masm");
  mesitylene.canonicalize();
  const auto mesitylenSubstitutionEdge = BondIndex {0, 14};

  Molecule nhc = IO::read("masm/nhc.masm");
  nhc.canonicalize();
  const auto nhcSubstitutionEdge = BondIndex {0, 7};

  Molecule substituted = Editing::substitute(
    nhc,
    mesitylene,
    nhcSubstitutionEdge,
    mesitylenSubstitutionEdge
  );

  BOOST_CHECK(substituted.graph().N() == mesitylene.graph().N() + nhc.graph().N() - 2);
}