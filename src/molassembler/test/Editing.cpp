/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include <boost/test/unit_test.hpp>
#include "boost/range/iterator_range_core.hpp"

#include "molassembler/Cycles.h"
#include "molassembler/Editing.h"
#include "molassembler/IO.h"
#include "molassembler/IO/SmilesParser.h"
#include "molassembler/Molecule.h"
#include "molassembler/Graph.h"
#include "molassembler/Subgraphs.h"

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
 * - multidentate_ligand: "CC(C)(C)P(CC1=NC(=CC=C1)CP(C(C)(C)C)C(C)(C)C)C(C)(C)C"
 * - haptic_ligand: "CC(C)(C)N[Si](C)(C)C1=C[C+]C=C1"
 *
 * The .cbor are just io::LineNotation::fromCanonicalSMILES and directly
 * exported using io::write(str, mol).
 */

using namespace Scine;
using namespace molassembler;

namespace detail {

auto makeIsBridgePredicate(const Molecule& mol) {
  return [&mol](const BondIndex& bond) -> bool {
    return !mol.graph().canRemove(bond);
  };
}

auto makeBondElementsPredicate(const Molecule& mol, const Utils::ElementType e1, const Utils::ElementType e2) {
  return [&mol, e1, e2](const BondIndex& bond) -> bool {
    const Graph& g = mol.graph();
    return (
      (g.elementType(bond.first) == e1 && g.elementType(bond.second) == e2)
      || (g.elementType(bond.first) == e2 && g.elementType(bond.second) == e1)
    );
  };
}

template<typename T, typename P1, typename P2>
auto combineUnaryPredicates(P1&& p1, P2&& p2) {
  return [&p1, &p2](const T& t) -> bool {
    return p1(t) && p2(t);
  };
}

boost::optional<AtomIndex> findSingle(const Molecule& mol, const Utils::ElementType element) {
  for(const AtomIndex i : boost::make_iterator_range(mol.graph().atoms())) {
    if(mol.graph().elementType(i) == element) {
      return {i};
    }
  }

  return boost::none;
}

std::vector<AtomIndex> findMultiple(const Molecule& mol, const Utils::ElementType element) {
  std::vector<AtomIndex> matches;

  for(const AtomIndex i : boost::make_iterator_range(mol.graph().atoms())) {
    if(mol.graph().elementType(i) == element) {
      matches.push_back(i);
    }
  }

  return matches;
}

template<typename UnaryPredicate>
boost::optional<BondIndex> findEdge(const Molecule& mol, UnaryPredicate&& predicate) {
  for(const BondIndex bond : boost::make_iterator_range(mol.graph().bonds())) {
    if(predicate(bond)) {
      return bond;
    }
  }

  return boost::none;
}

} // namespace detail

BOOST_AUTO_TEST_CASE(EditingCleave) {
  auto makeNMe = []() -> std::pair<Molecule, BondIndex> {
    Molecule methyl = io::experimental::parseSmilesSingleMolecule("[CH3]");
    std::vector<AtomIndex> methylPlugAtoms(1, 0);

    const AtomIndex CIndex = methylPlugAtoms.front();
    const AtomIndex NIndex = methyl.addAtom(Utils::ElementType::N, CIndex);
    return {
      methyl,
      BondIndex {CIndex, NIndex}
    };
  };

  Molecule caffeine = io::read("cbor/caffeine.cbor");
  caffeine.canonicalize();

  // Find a N-Me bridge bond to cleave
  const auto pattern = makeNMe();
  const auto matches = subgraphs::maximum(
    pattern.first,
    caffeine
  );
  BOOST_REQUIRE_MESSAGE(
    matches.size() > 0,
    "No matches found for N-Me pattern in caffeine!"
  );

  // One of the relevant matched N-C edges must be a bridge edge
  bool foundBridgeEdge = false;
  BondIndex bridgeEdge;
  for(const auto& match : matches) {
    bridgeEdge = {
      match.left.at(pattern.second.first),
      match.left.at(pattern.second.second)
    };

    if(!caffeine.graph().canRemove(bridgeEdge)) {
      foundBridgeEdge = true;
      break;
    }
  }

  BOOST_REQUIRE_MESSAGE(
    foundBridgeEdge,
    "Could not find a bridge edge from among matches for N-Me pattern in caffeine!"
  );

  // Perform the cleave test
  std::pair<Molecule, Molecule> cleaved;
  BOOST_REQUIRE_NO_THROW(cleaved = Editing::cleave(caffeine, bridgeEdge));

  BOOST_CHECK_MESSAGE(
    cleaved.first.graph().N() + cleaved.second.graph().N() == caffeine.graph().N(),
    "Sum of number of vertices of cleaved molecules does not match that of the original molecule"
  );

  BOOST_CHECK_MESSAGE(
    cleaved.first.graph().B() + cleaved.second.graph().B() == caffeine.graph().B() - 1,
    "Sum of number of bonds of cleaved molecules does not match that of the original molecule minus one"
  );
}

BOOST_AUTO_TEST_CASE(EditingInsert) {
  // Set up the test and make sure the canonical indices still match expectations
  Molecule biphenyl = io::read("cbor/biphenyl.cbor");
  biphenyl.canonicalize();

  // Find the C-C bridge edge
  auto bridgeOption = detail::findEdge(
    biphenyl,
    detail::combineUnaryPredicates<BondIndex>(
      detail::makeIsBridgePredicate(biphenyl),
      detail::makeBondElementsPredicate(biphenyl, Utils::ElementType::C, Utils::ElementType::C)
    )
  );
  BOOST_REQUIRE_MESSAGE(bridgeOption, "Could not find bridge edge in biphenyl");

  auto sides = biphenyl.graph().splitAlongBridge(*bridgeOption);

  BOOST_REQUIRE_MESSAGE(
    sides.first.size() == sides.second.size(),
    "Selected bridge edge does not split biphenyl equally in two"
  );

  Molecule pyrimidine = io::read("cbor/pyrimidine.cbor");
  pyrimidine.canonicalize();

  auto pyrimidineNitrogens = detail::findMultiple(pyrimidine, Utils::ElementType::N);

  BOOST_REQUIRE_MESSAGE(
    pyrimidineNitrogens.size() >= 2,
    "Could not find two nitrogens in pyrimidine"
  );

  // Insert pyrimidine into the bridge bond of biphenyl
  Molecule inserted = Editing::insert(
    biphenyl,
    pyrimidine,
    *bridgeOption,
    pyrimidineNitrogens.front(),
    pyrimidineNitrogens.at(1)
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

BOOST_AUTO_TEST_CASE(EditingSuperpose) {
  Molecule pyridine = io::read("cbor/pyridine.cbor");
  pyridine.canonicalize();

  auto pyridineNitrogenOption = detail::findSingle(pyridine, Utils::ElementType::N);
  BOOST_REQUIRE_MESSAGE(pyridineNitrogenOption, "Could not find N in pyridine");

  Molecule methane = io::read("cbor/methane.cbor");
  methane.canonicalize();

  auto methaneHydrogenOption = detail::findSingle(methane, Utils::ElementType::H);
  BOOST_REQUIRE_MESSAGE(methaneHydrogenOption, "Could not find H in methane");

  Molecule superposition = Editing::superpose(
    pyridine,
    methane,
    *pyridineNitrogenOption,
    *methaneHydrogenOption
  );

  BOOST_CHECK(superposition.graph().N() == pyridine.graph().N() + methane.graph().N() - 1);
}

BOOST_AUTO_TEST_CASE(EditingSubstitute) {
  Molecule chlorobenzene = io::read("cbor/chlorobenzene.cbor");
  chlorobenzene.canonicalize();

  auto chlorobenzeneBridgeOption = detail::findEdge(
    chlorobenzene,
    detail::combineUnaryPredicates<BondIndex>(
      detail::makeIsBridgePredicate(chlorobenzene),
      detail::makeBondElementsPredicate(
        chlorobenzene,
        Utils::ElementType::C,
        Utils::ElementType::Cl
      )
    )
  );
  BOOST_REQUIRE_MESSAGE(chlorobenzeneBridgeOption, "Could not find C-Cl bridge edge in chlorobenzene");

  Molecule phenole = io::read("cbor/phenole.cbor");
  phenole.canonicalize();
  auto phenoleBridgeOption = detail::findEdge(
    phenole,
    detail::combineUnaryPredicates<BondIndex>(
      detail::makeIsBridgePredicate(phenole),
      detail::makeBondElementsPredicate(
        phenole,
        Utils::ElementType::C,
        Utils::ElementType::O
      )
    )
  );
  BOOST_REQUIRE_MESSAGE(phenoleBridgeOption, "Could not find C-O bridge edge in phenole");

  Molecule substituted = Editing::substitute(
    chlorobenzene,
    phenole,
    *chlorobenzeneBridgeOption,
    *phenoleBridgeOption
  );

  Molecule biphenyl = io::read("cbor/biphenyl.cbor");
  BOOST_CHECK_MESSAGE(
    substituted == biphenyl,
    "Result of substitution not recognized as biphenyl prior to canonicalization"
  );

  substituted.canonicalize();
  biphenyl.canonicalize();
  BOOST_CHECK_MESSAGE(
    substituted == biphenyl,
    "Result of substitution not recognized as biphenyl after canonicalization"
  );
}

BOOST_AUTO_TEST_CASE(EditingConnect) {
  Molecule pyridine = io::read("cbor/pyridine.cbor");
  pyridine.canonicalize();
  auto pyridineNitrogenOption = detail::findSingle(pyridine, Utils::ElementType::N);
  BOOST_REQUIRE_MESSAGE(pyridineNitrogenOption, "Could not find N in pyridine");

  Molecule connected = Editing::connect(
    pyridine,
    pyridine,
    *pyridineNitrogenOption,
    *pyridineNitrogenOption,
    BondType::Single
  );

  BOOST_CHECK(connected.graph().N() == 2 * pyridine.graph().N());
}

BOOST_AUTO_TEST_CASE(EditingBugfixMesityleneSubstitution) {
  Molecule mesitylene = io::read("cbor/mesitylene.cbor");
  mesitylene.canonicalize();

  const auto mesityleneSubstitutionEdgeOption = detail::findEdge(
    mesitylene,
    detail::combineUnaryPredicates<BondIndex>(
      detail::makeIsBridgePredicate(mesitylene),
      detail::makeBondElementsPredicate(
        mesitylene,
        Utils::ElementType::C,
        Utils::ElementType::H
      )
    )
  );
  BOOST_REQUIRE_MESSAGE(
    mesityleneSubstitutionEdgeOption,
    "Could not find C-C bridge edge in mesitylene"
  );

  Molecule nhc = io::read("cbor/nhc.cbor");
  nhc.canonicalize();

  Molecule pattern {
    Utils::ElementType::C,
    Utils::ElementType::H
  };
  for(unsigned i = 0; i < 2; ++i) {
    pattern.addAtom(Utils::ElementType::N, 0);
  }

  auto matches = subgraphs::maximum(pattern, nhc);
  BOOST_REQUIRE_MESSAGE(!matches.empty(), "Could not find C(HNN) pattern in NHC");

  // Map C-H bond indices from pattern to nhc
  const BondIndex nhcSubstitutionEdge {
    matches.front().left.at(0),
    matches.front().left.at(1)
  };

  Molecule substituted = Editing::substitute(
    nhc,
    mesitylene,
    nhcSubstitutionEdge,
    *mesityleneSubstitutionEdgeOption
  );

  BOOST_CHECK(substituted.graph().N() == mesitylene.graph().N() + nhc.graph().N() - 2);
}

BOOST_AUTO_TEST_CASE(EditingBugfixHapticLigands) {
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

BOOST_AUTO_TEST_CASE(EditingAddMultidentateLigand) {
  Molecule ligand = io::read("cbor/multidentate_ligand.cbor");
  ligand.canonicalize();

  auto NOption = detail::findSingle(ligand, Utils::ElementType::N);
  auto Ps = detail::findMultiple(ligand, Utils::ElementType::P);

  BOOST_REQUIRE(NOption && Ps.size() == 2);

  std::vector<AtomIndex> ligandBindingAtoms = Ps;
  ligandBindingAtoms.push_back(*NOption);

  Molecule complex {
    Utils::ElementType::Ru,
    Utils::ElementType::H
  };

  BOOST_CHECK_NO_THROW(
    complex = Editing::addLigand(complex, ligand, 0, ligandBindingAtoms)
  );
}

BOOST_AUTO_TEST_CASE(EditingAddHapticLigand) {
  auto ligand = io::read("cbor/haptic_ligand.cbor");

  std::vector<AtomIndex> ligandCycle;
  for(const auto& cycle : ligand.graph().cycles()) {
    ligandCycle = makeRingIndexSequence(cycle);
  }

  BOOST_REQUIRE_MESSAGE(ligandCycle.size() == 5, "Could not find cycle in haptic_ligand.cbor");

  auto NOption = detail::findSingle(ligand, Utils::ElementType::N);
  BOOST_REQUIRE_MESSAGE(NOption, "Could not find nitrogen in haptic_ligand.cbor");

  auto ligandBindingAtoms = ligandCycle;
  ligandBindingAtoms.push_back(*NOption);

  Molecule mol {Utils::ElementType::Ti};
  for(unsigned i = 0; i < 2; ++i) {
    mol.addAtom(Utils::ElementType::Cl, 0);
  }

  Molecule complex;
  BOOST_CHECK_NO_THROW(
    complex = Editing::addLigand(mol, ligand, 0, ligandBindingAtoms)
  );
}
