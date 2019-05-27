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
 * The .cbor are just IO::LineNotation::fromCanonicalSMILES and directly
 * exported using IO::write(str, mol).
 */

BOOST_AUTO_TEST_CASE(cleaveEdit) {
  using namespace Scine;
  using namespace molassembler;

  Molecule caffeine = IO::read("cbor/caffeine.cbor");
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
  Molecule biphenyl = IO::read("cbor/biphenyl.cbor");
  biphenyl.canonicalize();

  const BondIndex bridge {12, 13};
  BOOST_REQUIRE_MESSAGE(!biphenyl.graph().canRemove(bridge), "Selected edge for test is not a bridge edge");
  auto sides = biphenyl.graph().splitAlongBridge(bridge);

  BOOST_REQUIRE_MESSAGE(
    sides.first.size() == sides.second.size(),
    "Selected bridge edge does not split biphenyl equally in two"
  );

  Molecule pyrimidine = IO::read("cbor/pyrimidine.cbor");
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

  Molecule pyridine = IO::read("cbor/pyridine.cbor");
  pyridine.canonicalize();
  const AtomIndex pyridineNitrogen = 5;
  BOOST_REQUIRE(pyridine.graph().elementType(pyridineNitrogen) == Utils::ElementType::N);

  Molecule methane = IO::read("cbor/methane.cbor");
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

  Molecule chlorobenzene = IO::read("cbor/chlorobenzene.cbor");
  chlorobenzene.canonicalize();
  const BondIndex chlorobenzeneBridge {5, 7};
  BOOST_REQUIRE_MESSAGE(!chlorobenzene.graph().canRemove(chlorobenzeneBridge), "Selected edge for test is not a bridge edge");

  Molecule phenole = IO::read("cbor/phenole.cbor");
  phenole.canonicalize();
  const BondIndex phenoleBridge {6, 8};
  BOOST_REQUIRE_MESSAGE(!phenole.graph().canRemove(phenoleBridge), "Selected edge for test is not a bridge edge");

  Molecule substituted = Editing::substitute(
    chlorobenzene,
    phenole,
    chlorobenzeneBridge,
    phenoleBridge
  );

  Molecule biphenyl = IO::read("cbor/biphenyl.cbor");
  BOOST_CHECK(substituted == biphenyl);

  substituted.canonicalize();
  biphenyl.canonicalize();
  BOOST_CHECK(substituted == biphenyl);
}

BOOST_AUTO_TEST_CASE(connectEdit) {
  using namespace Scine;
  using namespace molassembler;

  Molecule pyridine = IO::read("cbor/pyridine.cbor");
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

  Molecule mesitylene = IO::read("cbor/mesitylene.cbor");
  mesitylene.canonicalize();
  const auto mesitylenSubstitutionEdge = BondIndex {0, 14};

  Molecule nhc = IO::read("cbor/nhc.cbor");
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

BOOST_AUTO_TEST_CASE(makingHapticLigandTest) {
  using namespace Scine;
  using namespace molassembler;

  Molecule complex {
    Utils::ElementType::Ru,
    Utils::ElementType::H
  };

  { // Variant one: Extend a hydrogen atom to hydrogen molecule, then complexate
    Molecule complexCopy = complex;
    AtomIndex hydrogenOne = complexCopy.addAtom(Utils::ElementType::H, 0, BondType::Single);
    AtomIndex hydrogenTwo = complexCopy.addAtom(Utils::ElementType::H, hydrogenOne, BondType::Single);
    BOOST_CHECK_NO_THROW(complexCopy.addBond(0, hydrogenTwo, BondType::Single));
  }

  { // Variant two: Connect two individual bonding hydrogen atoms
    /*Molecule complexCopy = complex;
    AtomIndex hydrogenOne = complexCopy.addAtom(Utils::ElementType::H, 0, BondType::Single);
    AtomIndex hydrogenTwo = complexCopy.addAtom(Utils::ElementType::H, 0, BondType::Single);
    complexCopy.addBond(hydrogenOne, hydrogenTwo, BondType::Single);*/
  }
}

BOOST_AUTO_TEST_CASE(addLigandTest) {
  using namespace Scine;
  using namespace molassembler;

  Molecule ligand = IO::LineNotation::fromCanonicalSMILES("CC(C)(C)P(CC1=NC(=CC=C1)CP(C(C)(C)C)C(C)(C)C)C(C)(C)C");
  ligand.canonicalize();
  BOOST_CHECK(ligand.graph().elementType(43) == Utils::ElementType::N);
  BOOST_CHECK(ligand.graph().elementType(49) == Utils::ElementType::P);
  BOOST_CHECK(ligand.graph().elementType(50) == Utils::ElementType::P);

  Molecule complex {
    Utils::ElementType::Ru,
    Utils::ElementType::H
  };

  BOOST_CHECK_NO_THROW(
    complex = Editing::addLigand(
      complex,
      ligand,
      0,
      {43, 49, 50}
    )
  );
}

BOOST_AUTO_TEST_CASE(hapticLigandAddition) {
  using namespace Scine;
  using namespace molassembler;

  auto ligand = IO::LineNotation::fromCanonicalSMILES("CC(C)(C)N[Si](C)(C)C1=C[C+]C=C1");
  Molecule mol {Utils::ElementType::Ti};
  for(unsigned i = 0; i < 2; ++i) {
    mol.addAtom(Utils::ElementType::Cl, 0);
  }

  Molecule complex;
  BOOST_CHECK_NO_THROW(
    complex = Editing::addLigand(
      mol,
      ligand,
      0,
      {4, 8, 9, 10, 11, 12}
    )
  );
}
