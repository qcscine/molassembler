/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_TEST_MODULE MolassemblerMainTests
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/test/unit_test.hpp"

#include "shapes/Data.h"

#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"
#include "temple/Random.h"
#include "temple/Invoke.h"
#include "temple/Stringify.h"
#include "temple/Optionals.h"

#include "molassembler/Graph.h"
#include "molassembler/Graph/PrivateGraph.h"
#include "molassembler/IO.h"
#include "molassembler/Interpret.h"
#include "molassembler/Isomers.h"
#include "molassembler/Molecule.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/Options.h"
#include "molassembler/AtomStereopermutator.h"
#include "molassembler/BondStereopermutator.h"
#include "molassembler/StereopermutatorList.h"

#include "Utils/Geometry/ElementInfo.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/IO/ChemicalFileFormats/ChemicalFileHandler.h"

#include "temple/UnorderedSetAlgorithms.h"
#include <iostream>

using namespace Scine;
using namespace molassembler;

static_assert(
  AtomEnvironmentComponents::All == (
    AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Shapes
    | AtomEnvironmentComponents::Stereopermutations
  ),
  "All needs to include all components in its definition"
);

std::tuple<Molecule, Molecule, std::vector<AtomIndex>> readIsomorphism(const boost::filesystem::path& filePath) {
  auto readData = Utils::ChemicalFileHandler::read(filePath.string());
  auto permutedData = io::shuffle(readData.first, readData.second);

  auto interpretSingle = [](const Utils::AtomCollection& ac, const Utils::BondOrderCollection& boc) -> Molecule {
    interpret::MoleculesResult interpretation;
    if(boc.empty()) {
      // Unfortunately, the file type does not include bond order information
      interpretation = interpret::molecules(ac, interpret::BondDiscretizationOption::RoundToNearest);
    } else {
      interpretation = interpret::molecules(ac, boc, interpret::BondDiscretizationOption::RoundToNearest);
    }

    if(interpretation.molecules.size() > 1) {
      throw std::runtime_error(
        std::string("File is not a single molecule, but contains ")
          + std::to_string(interpretation.molecules.size())
          + " components."
      );
    }

    return interpretation.molecules.front();
  };


  return std::make_tuple(
    interpretSingle(readData.first, readData.second),
    interpretSingle(std::get<0>(permutedData), std::get<1>(permutedData)),
    std::move(std::get<2>(permutedData))
  );
}

// Molecule instances are trivial to handle and can do all they should
BOOST_AUTO_TEST_CASE(MoleculeRuleOfFiveTrivial) {
  // Default constructor
  Molecule f, g;

  // Copy assignment
  f = g;

  // Move assignment
  g = Molecule {};

  // Copy constructor
  Molecule h {g};

  // Move constructor
  Molecule i {Molecule {}};

  std::vector<Molecule> allJustHydrogen {f, g, h, i};

  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::adaptors::allPairs(allJustHydrogen),
      std::equal_to<>()
    ),
    "Basic construction does not generate equal molecules"
  );
}

using HashArgumentsType = std::tuple<
  Utils::ElementType,
  std::vector<molassembler::hashes::BondInformation>,
  boost::optional<shapes::Shape>,
  boost::optional<unsigned>
>;

std::string repr(const HashArgumentsType& hashArgs) {
  std::string representation;
  representation += Utils::ElementInfo::symbol(std::get<0>(hashArgs));
  representation += ", ";
  representation += temple::condense(
    temple::map(
      std::get<1>(hashArgs),
      [](const hashes::BondInformation& b) -> std::string {
        return (
          "bty " + std::to_string(static_cast<unsigned>(b.bondType))
          + ", s = " + std::to_string(b.stereopermutatorOnBond)
          + ", a = " + temple::optionals::map(
            b.assignmentOptional,
            [](const unsigned i) -> std::string {
              return std::to_string(i);
            }
          ).value_or("N")
        );
      }
    )
  );
  representation += ", s = ";
  representation += temple::optionals::map(
    std::get<2>(hashArgs),
    [](const shapes::Shape s) -> std::string {
      return shapes::name(s);
    }
  ).value_or("N");
  representation += ", a = ";
  representation += temple::optionals::map(
    std::get<3>(hashArgs),
    [](const unsigned p) -> std::string {
      return std::to_string(p);
    }
  ).value_or("N");
  return representation;
}

HashArgumentsType randomArguments() {
  auto genBondInformation = []() -> molassembler::hashes::BondInformation {
    auto bty = static_cast<BondType>(
      temple::random::getSingle<unsigned>(0, 6, randomnessEngine())
    );

    bool bondStereopermutatorPresent = temple::random::getSingle<bool>(randomnessEngine());
    boost::optional<unsigned> bondStereopermutatorAssignment = boost::none;

    if(bondStereopermutatorPresent) {
      bondStereopermutatorAssignment = temple::random::getSingle<unsigned>(0, 1, randomnessEngine());
    }

    return {
      bty,
      bondStereopermutatorPresent,
      bondStereopermutatorAssignment
    };
  };

  boost::optional<shapes::Shape> shapeOptional;
  boost::optional<unsigned> assignmentOptional;
  if(temple::random::getSingle<bool>(randomnessEngine())) {
    shapeOptional = static_cast<shapes::Shape>(
      temple::random::getSingle<unsigned>(0, 15, randomnessEngine())
    );

    std::geometric_distribution<unsigned> gd {0.2};
    assignmentOptional = gd(randomnessEngine());
  }

  // If a symmetry is specified, the bond number must match
  std::vector<molassembler::hashes::BondInformation> bonds;
  unsigned S;
  if(shapeOptional) {
    S = shapes::size(*shapeOptional);
  } else {
    S = temple::random::getSingle<unsigned>(1, 8, randomnessEngine());
  }

  for(unsigned i = 0; i < S; ++i) {
    bonds.emplace_back(genBondInformation());
  }
  std::sort(bonds.begin(), bonds.end());

  return {
    Utils::ElementInfo::element(
      temple::random::getSingle<unsigned>(1, 112, randomnessEngine())
    ),
    bonds,
    shapeOptional,
    assignmentOptional
  };
}

// Atom environment wide hashes cannot be collided
BOOST_AUTO_TEST_CASE(AtomEnvironmentHashesDoNotCollide) {
  auto bitmaskTuple = std::make_tuple(
    AtomEnvironmentComponents::All
  );

  // Try to guess a disjoint combination that has the same value
  std::unordered_map<molassembler::hashes::WideHashType, HashArgumentsType> resultsMap;
  for(unsigned N = 0; N < 1e6; ++N) {
    auto arguments = randomArguments();

    auto result = temple::detail::invokeHelper(
      hashes::hash,
      std::tuple_cat(bitmaskTuple, arguments),
      std::make_index_sequence<5> {}
    );

    auto findIter = resultsMap.find(result);
    if(findIter != std::end(resultsMap)) {
      BOOST_REQUIRE_MESSAGE(
        arguments == findIter->second,
        "Found overlapping result for different arguments to hashAtomEnvironment!"
        << "\nA: " << repr(arguments) << "\nB: " << repr(findIter->second)
      );
    } else {
      resultsMap.emplace(result, std::move(arguments));
    }
  }
}

/* Hashes, stereopermutator lists and graphs are identical across two molecules
 * generated by applying a random permutation to the intermediate data (atoms
 * and BOs) and applying the same permutation via .applyPermutation(perm)
 */
BOOST_AUTO_TEST_CASE(MoleculeGraphPermutation) {
  boost::filesystem::path directoryBase("isomorphisms");

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    if(currentFilePath.extension() != ".mol") {
      continue;
    }

    Molecule a, b;
    std::vector<AtomIndex> permutation;
    std::tie(a, b, permutation) = readIsomorphism(currentFilePath);

    Molecule permuted = a;
    permuted.applyPermutation(permutation);

    // b and permuted must be 1:1 identical
    auto bHashes = hashes::generate(b.graph().inner(), b.stereopermutators(), AtomEnvironmentComponents::All);
    auto permutedHashes = hashes::generate(permuted.graph().inner(), permuted.stereopermutators(), AtomEnvironmentComponents::All);

    BOOST_CHECK_MESSAGE(
      bHashes == permutedHashes,
      "Hashes for " << currentFilePath.string() << " permuted in two fashions do not match"
    );

    BOOST_CHECK_MESSAGE(
      b.stereopermutators() == permuted.stereopermutators(),
      "Stereopermutator lists for " << currentFilePath.string() << " permuted in two fashions do not match"
    );

    BOOST_CHECK_MESSAGE(
      b.graph().inner().identicalGraph(permuted.graph().inner()),
      "Graphs for " << currentFilePath.string() << " permuted in two fashions do not match"
    );
  }
}

/* Hashes are identical after applying a permutation (this is also tested in
 * moleculePermutation in the inverse fashion, I think, and could be removed)
 */
BOOST_AUTO_TEST_CASE(AtomEnvironmentHashesRegular) {
  boost::filesystem::path directoryBase("isomorphisms");

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    if(currentFilePath.extension() != ".mol") {
      continue;
    }

    Molecule a = io::read(currentFilePath.string());

    auto aWideHashes = hashes::generate(a.graph().inner(), a.stereopermutators(), AtomEnvironmentComponents::All);

    auto permutation = temple::iota<AtomIndex>(a.graph().N());
    temple::random::shuffle(permutation, randomnessEngine());

    Molecule b = a;
    b.applyPermutation(permutation);

    auto bWideHashes = hashes::generate(b.graph().inner(), b.stereopermutators(), AtomEnvironmentComponents::All);

    for(AtomIndex i = 0; i < a.graph().N(); ++i) {
      BOOST_CHECK_MESSAGE(
        aWideHashes.at(i) == bWideHashes.at(permutation.at(i)),
        "Mismatch: hash(a, " << i << ") = " << aWideHashes.at(i) << " != " << bWideHashes.at(permutation.at(i)) << " = hash(b, " << permutation.at(i) << ")"
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(MoleculeCanonicalizationAtomMap) {
  boost::filesystem::path directoryBase("ranking_tree_molecules");

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    if(currentFilePath.extension() != ".mol") {
      continue;
    }

    const Molecule m = io::read(currentFilePath.string());
    Molecule n = m;
    auto indexMap = n.canonicalize();

    BOOST_CHECK_MESSAGE(
      temple::all_of(
        temple::adaptors::range(m.graph().N()),
        [&](const unsigned oldIndex) -> bool {
          unsigned newIndex = indexMap.at(oldIndex);
          return n.graph().elementType(newIndex) == m.graph().elementType(oldIndex);
        }
      ),
      "Index map from canonicalization is not consistent with reordering"
    );
  }
}

// After canonicalizing isomorphic molecules, they are identical
BOOST_AUTO_TEST_CASE(MoleculeCanonicalization) {
  boost::filesystem::path directoryBase("isomorphisms");

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    if(currentFilePath.extension() != ".mol") {
      continue;
    }

    Molecule a, b;
    std::tie(a, b, std::ignore) = readIsomorphism(currentFilePath);

    // Canonicalize both molecules
    a.canonicalize();
    b.canonicalize();

    // These must be IDENTICAL afterwards, not isomorphic
    auto aWideHashes = hashes::generate(a.graph().inner(), a.stereopermutators(), AtomEnvironmentComponents::All);
    auto bWideHashes = hashes::generate(b.graph().inner(), b.stereopermutators(), AtomEnvironmentComponents::All);

    BOOST_CHECK_MESSAGE(
      aWideHashes == bWideHashes,
      "After canonizing two isomorphic instances of " << currentFilePath.stem()
      << ", their hashes are not lexicographically equal"
    );

    BOOST_CHECK_MESSAGE(
      a.graph().inner().identicalGraph(b.graph().inner()),
      "After canonizing two isomorphic instances of " << currentFilePath.stem()
      << ", their graphs are not identical (1:1 the same, not an isomorphism test)"
    );
  }
}

BOOST_AUTO_TEST_CASE(MoleculeHashes) {
  boost::filesystem::path directoryBase("isomorphisms");

  using C = AtomEnvironmentComponents;
  std::vector<C> testComponents {
    C::Connectivity,
    C::Connectivity | C::BondOrders,
    C::Connectivity | C::BondOrders | C::Shapes,
    C::All
  };

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    if(currentFilePath.extension() != ".mol") {
      continue;
    }

    Molecule a, b;
    std::tie(a, b, std::ignore) = readIsomorphism(currentFilePath);

    for(C components : testComponents) {
      Molecule c = a;
      c.canonicalize(components);
      Molecule d = b;
      d.canonicalize(components);

      BOOST_CHECK_EQUAL(c.hash(), d.hash());
    }
  }
}

// Isomorphic molecules are recognized as such by modularCompare
BOOST_AUTO_TEST_CASE(MoleculeIsomorphism) {
  using namespace std::string_literals;

  boost::filesystem::path directoryBase("isomorphisms");

  std::vector<Molecule> originals;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    if(currentFilePath.extension() != ".mol") {
      continue;
    }

    Molecule a, b;
    std::tie(a, b, std::ignore) = readIsomorphism(currentFilePath);

    bool pass = (a == b);

    BOOST_CHECK_MESSAGE(
      pass,
      "Molecule isomorphism fails for " << currentFilePath.string() << "!"
    );

    if(!pass) {
      boost::filesystem::path base = currentFilePath.stem();
      std::string aPathString = base.string() + ".dot";
      std::string bPathString = base.string() + "isomorphic_.dot";
      std::cout << "Writing dot files '" << aPathString << "'.\n";
      std::ofstream aGraph(aPathString);
      aGraph << a.dumpGraphviz();
      aGraph.close();
      std::ofstream bGraph(bPathString);
      bGraph << b.dumpGraphviz();
      aGraph.close();
    }

    originals.push_back(std::move(a));
  }

  BOOST_REQUIRE_MESSAGE(
    !originals.empty(),
    "No molecules found for isomorphism checks"
  );

  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::adaptors::allPairs(originals),
      std::not_equal_to<>()
    ),
    "Some originals in the isomorphism test folder match one another!"
  );
}

BOOST_AUTO_TEST_CASE(MoleculeIsotopicInequivalency) {
  Molecule H2;
  Molecule HD;
  HD.setElementType(0, Utils::ElementType::D);
  BOOST_CHECK(H2 != HD);
}

// Atom stereocenter assignments are part of strict equivalency
BOOST_AUTO_TEST_CASE(MoleculeBasicRSInequivalency) {
  // Build an asymmetric tetrahedral carbon
  Molecule a {Utils::ElementType::C, Utils::ElementType::H, BondType::Single};
  a.addAtom(Utils::ElementType::F, 0, BondType::Single);
  a.addAtom(Utils::ElementType::Cl, 0, BondType::Single);
  a.addAtom(Utils::ElementType::Br, 0, BondType::Single);
  a.setShapeAtAtom(0, shapes::Shape::Tetrahedron);

  // Make sure it's recognized as asymmetric
  auto centralStereopermutatorOption = a.stereopermutators().option(0);
  BOOST_CHECK(
    centralStereopermutatorOption
    && centralStereopermutatorOption->numStereopermutations() == 2
  );

  // Assign it and create its opposite stereopermutation in another Molecule
  a.assignStereopermutator(0, 0);
  Molecule b = a;
  b.assignStereopermutator(0, 1);

  // These must compare unequal
  BOOST_CHECK(a != b);
}

// Bond stereocenter assignments are part of strict equivalency
BOOST_AUTO_TEST_CASE(MoleculeBasicEZInequivalency) {
  Molecule a {Utils::ElementType::C, Utils::ElementType::C, BondType::Double};
  a.addAtom(Utils::ElementType::H, 0, BondType::Single);
  a.addAtom(Utils::ElementType::F, 0, BondType::Single);
  a.addAtom(Utils::ElementType::H, 1, BondType::Single);
  a.addAtom(Utils::ElementType::F, 1, BondType::Single);

  // Set the shapes
  a.setShapeAtAtom(0, shapes::Shape::EquilateralTriangle);
  a.setShapeAtAtom(1, shapes::Shape::EquilateralTriangle);

  // Progression must recognize the new stereopermutator
  auto stereopermutatorOption = a.stereopermutators().option(
    BondIndex {0, 1}
  );
  BOOST_CHECK(
    stereopermutatorOption
    && stereopermutatorOption->numStereopermutations() == 2
  );

  a.assignStereopermutator(
    BondIndex {0, 1},
    0
  );

  Molecule b = a;
  b.assignStereopermutator(
    BondIndex {0, 1},
    1
  );

  // These must compare unequal
  BOOST_CHECK(a != b);
}

// Supplied by RankingSpotCheck
extern bool isStereogenic(
  const Molecule& molecule,
  AtomIndex i
);

/* If a stereocenter is changed, ranking effects are propagated correctly
 * at other stereocenters.
 */
BOOST_AUTO_TEST_CASE(PropagateGraphChangeTests) {
  boost::filesystem::path filePath("ranking_tree_molecules");
  filePath /= "(2R,3r,4S)-pentane-2,3,4-trithiol.mol";

  auto pseudocenter = io::read(
    filePath.string()
  );

  AtomIndex central = 0;
  std::array<AtomIndex, 2> outer {{1, 6}};

  /* If the outer stereopermutators have the same assignment, the central
   * stereopermutator shouldn't exist. If the outer stereopermutators have a
   * different assignment, the central stereopermutator has to exist.
   *
   * Initially, we have
   * 0 -> 1 r
   * 1 -> 1 R
   * 6 -> 0 S
   */
  // Make 1 from R to S -> stereopermutator should disappear
  pseudocenter.assignStereopermutator(outer.front(), 0);

  BOOST_CHECK(!isStereogenic(pseudocenter, central));

  // Make 6 from S to R -> stereopermutator should reappear
  pseudocenter.assignStereopermutator(outer.back(), 1);

  BOOST_CHECK(isStereogenic(pseudocenter, central));


  /* Bug reported by Stephanie: Central stereopermutator does not propagate
   * to original abstract site case after forcing cut-tetrahedral shapes
   * for each nitrogen.
   */
  auto complex = io::read("various/propagation-test-case-1.json");

  for(AtomIndex i = 0; i < complex.graph().N(); ++i) {
    if(
      complex.graph().elementType(i) == Utils::ElementType::N
      && complex.stereopermutators().option(i)
      && shapes::size(complex.stereopermutators().option(i)->getShape()) == 3
    ) {
      complex.setShapeAtAtom(i, shapes::Shape::VacantTetrahedron);
    }
  }

  auto stereopermutatorOption = complex.stereopermutators().option(0);
  BOOST_REQUIRE_MESSAGE(
    stereopermutatorOption,
    "The stereopermutator is no longer present!"
  );
  BOOST_CHECK_MESSAGE(
    stereopermutatorOption->numAssignments() > 2,
    "The stereopermutator does not have more than two feasible assignments, it has only " << stereopermutatorOption->numAssignments()
  );
}

BOOST_AUTO_TEST_CASE(MoleculeSplitRecognition) {
  std::vector<Molecule> molSplat, xyzSplat;

  BOOST_REQUIRE_NO_THROW(molSplat = io::split("multiple_molecules/multi_interpret.mol"));
  BOOST_REQUIRE_NO_THROW(xyzSplat = io::split("multiple_molecules/multi_interpret.xyz"));

  // Each file contains one ethane, four water, and a potassium

  Molecule water {Utils::ElementType::O};
  water.addAtom(Utils::ElementType::H, 0, BondType::Single);
  water.addAtom(Utils::ElementType::H, 0, BondType::Single);
  water.setShapeAtAtom(0, shapes::Shape::Bent);

  Molecule ethane {Utils::ElementType::C};
  ethane.addAtom(Utils::ElementType::H, 0, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, 0, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, 0, BondType::Single);
  AtomIndex otherCarbon = ethane.addAtom(Utils::ElementType::C, 0, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, otherCarbon, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, otherCarbon, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, otherCarbon, BondType::Single);
  ethane.setShapeAtAtom(0, shapes::Shape::Tetrahedron);
  ethane.setShapeAtAtom(otherCarbon, shapes::Shape::Tetrahedron);

  Molecule potassium {Utils::ElementType::K};

  auto countMolecule = [&](const std::vector<Molecule>& molecules, const Molecule& countMol) -> unsigned {
    return temple::accumulate(
      molecules,
      0u,
      [&](const unsigned carry, const Molecule& testMol) -> unsigned {
        if(testMol == countMol) {
          return carry + 1;
        }

        return carry;
      }
    );
  };

  BOOST_CHECK_MESSAGE(molSplat.size() == 6, "Mol separation yielded not 6, but " << molSplat.size() << " molecules");
  BOOST_CHECK_MESSAGE(xyzSplat.size() == 6, "XYZ separation yielded not 6, but " << xyzSplat.size() << " molecules");

  auto checkSplat = [&](const std::vector<Molecule>& splat) -> bool {
    unsigned waters = countMolecule(splat, water);
    unsigned ethanes = countMolecule(splat, ethane);
    unsigned potassiums = countMolecule(splat, potassium);

    BOOST_CHECK_MESSAGE(waters == 4u, "Did not get expected number of waters");
    BOOST_CHECK_MESSAGE(ethanes == 1u, "Did not get expected number of ethanes");
    BOOST_CHECK_MESSAGE(potassiums == 1u, "Did not get expected number of potassiums");

    return (waters == 4u && ethanes == 1u && potassiums == 1u);
  };

  BOOST_CHECK_MESSAGE(checkSplat(molSplat), "Molecule counts for interpret of multi_interpret.mol failed");
  BOOST_CHECK_MESSAGE(checkSplat(xyzSplat), "Molecule counts for interpret of multi_interpret.xyz failed");
}

BOOST_AUTO_TEST_CASE(MoleculeGeometryChoices) {
  molassembler::Molecule testMol(Utils::ElementType::Ru, Utils::ElementType::N, BondType::Single);
  testMol.addAtom(Utils::ElementType::H, 1u, BondType::Single);
  testMol.addAtom(Utils::ElementType::H, 1u, BondType::Single);


  auto stereocenterOption = testMol.stereopermutators().option(1u);
  BOOST_REQUIRE(stereocenterOption);

  if(auto suggestedShapeOption = testMol.inferShape(1u, stereocenterOption->getRanking())) {
    BOOST_CHECK(suggestedShapeOption.value() == shapes::Shape::VacantTetrahedron);
    testMol.setShapeAtAtom(1u, suggestedShapeOption.value());
  }

  testMol.addAtom(Utils::ElementType::H, 1u, BondType::Single);

  BOOST_CHECK(
    testMol.stereopermutators().option(1u)
    && testMol.stereopermutators().option(1u)->getShape() == shapes::Shape::Tetrahedron
  );
}

// Isomer predicates work as expected
BOOST_AUTO_TEST_CASE(IsomerPredicateTests) {
  Molecule a, b;
  BOOST_REQUIRE_NO_THROW(a = molassembler::io::read("isomers/enantiomers/Citalopram-R.mol"));
  BOOST_REQUIRE_NO_THROW(b = molassembler::io::read("isomers/enantiomers/Citalopram-S.mol"));

  constexpr auto bitmask = AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Shapes;

  a.canonicalize(bitmask);
  b.canonicalize(bitmask);

  BOOST_CHECK_MESSAGE(
    molassembler::enantiomeric(a, b),
    "Citalopram-R and -S are falsely determined not to be enantiomers."
  );

  BOOST_REQUIRE_NO_THROW(a = molassembler::io::read("isomers/enantiomers/Isoleucine-RS.mol"));
  BOOST_REQUIRE_NO_THROW(b = molassembler::io::read("isomers/enantiomers/Isoleucine-SR.mol"));

  a.canonicalize(bitmask);
  b.canonicalize(bitmask);

  BOOST_CHECK_MESSAGE(
    molassembler::enantiomeric(a, b),
    "Isoleucine-RS and -SR are falsely determined not to be enantiomers."
  );
}

BOOST_AUTO_TEST_CASE(EtaBondDynamism) {
  Molecule mol {Utils::ElementType::Fe, Utils::ElementType::C, BondType::Single};

  const AtomIndex carbonSubstituent = 1u;
  const AtomIndex secondCarbon = mol.addAtom(Utils::ElementType::C, carbonSubstituent, BondType::Double);
  mol.addAtom(Utils::ElementType::H, carbonSubstituent, BondType::Single);
  mol.addAtom(Utils::ElementType::H, carbonSubstituent, BondType::Single);
  mol.addAtom(Utils::ElementType::H, secondCarbon, BondType::Single);
  mol.addAtom(Utils::ElementType::H, secondCarbon, BondType::Single);

  auto ironSecondCarbonBond = mol.addBond(0, secondCarbon, BondType::Single);

  BOOST_CHECK(mol.graph().bondType(ironSecondCarbonBond) == BondType::Eta);
  BOOST_CHECK(mol.graph().bondType(BondIndex {0, 1}) == BondType::Eta);

  mol.removeBond(ironSecondCarbonBond);

  BOOST_CHECK(mol.graph().bondType(BondIndex {0, 1}) == BondType::Single);
}

void checkAtomStereopermutator(
  const Molecule& m,
  const AtomIndex i,
  const shapes::Shape shape
) {
  BOOST_CHECK_MESSAGE(
    temple::optionals::map(
      m.stereopermutators().option(i),
      [&](const AtomStereopermutator& permutator) {
        return permutator.getShape() == shape;
      }
    ).value_or(false),
    temple::optionals::map(
      m.stereopermutators().option(i),
      [&](const AtomStereopermutator& permutator) {
        return shapes::name(permutator.getShape());
      }
    ).value_or("No") << " atom stereopermutator on " << i << ", expected a(n) " << shapes::name(shape) << " stereopermutator"
  );
}

BOOST_AUTO_TEST_CASE(ShapeClassification) {
  // Particular cases from issues that reveal shape classification needs.

  // This is a huge system with octahedral iron centers (#90)
  auto interconnected = io::read("shape_classification/interconnected.mol");
  for(const AtomIndex i : interconnected.graph().atoms()) {
    if(interconnected.graph().elementType(i) == Utils::ElementType::Fe) {
      checkAtomStereopermutator(interconnected, i, shapes::Shape::Octahedron);
    }
  }

  // Trigonal bipyramidal transition metal center (#99)
  auto trigbipy = io::read("shape_classification/trig_bipy.mol");
  checkAtomStereopermutator(trigbipy, 0, shapes::Shape::TrigonalBipyramid);
}
