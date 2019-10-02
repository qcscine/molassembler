/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_TEST_MODULE MolassemblerMainTests
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/test/unit_test.hpp"

#include "chemical_symmetries/Symmetries.h"

#include "temple/Adaptors/AllPairs.h"
#include "temple/Adaptors/Iota.h"
#include "temple/Functional.h"
#include "temple/Random.h"
#include "temple/Invoke.h"
#include "temple/Stringify.h"

#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/IO.h"
#include "molassembler/Interpret.h"
#include "molassembler/Isomers.h"
#include "molassembler/Molecule.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/Options.h"
#include "molassembler/StereopermutatorList.h"

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
    | AtomEnvironmentComponents::Symmetries
    | AtomEnvironmentComponents::Stereopermutations
  ),
  "All needs to include all components in its definition"
);

std::tuple<Molecule, Molecule, std::vector<AtomIndex>> readIsomorphism(const boost::filesystem::path& filePath) {
  auto readData = Utils::ChemicalFileHandler::read(filePath.string());
  auto permutedData = IO::shuffle(readData.first, readData.second);

  auto interpretSingle = [](const Utils::AtomCollection& ac, const Utils::BondOrderCollection& boc) -> Molecule {
    InterpretResult interpretation;
    if(boc.empty()) {
      // Unfortunately, the file type does not include bond order information
      interpretation = interpret(ac, BondDiscretizationOption::RoundToNearest);
    } else {
      interpretation = interpret(ac, boc, BondDiscretizationOption::RoundToNearest);
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
  boost::optional<Symmetry::Name>,
  boost::optional<unsigned>
>;

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

  boost::optional<Symmetry::Name> symmetryOptional;
  boost::optional<unsigned> assignmentOptional;
  if(temple::random::getSingle<bool>(randomnessEngine())) {
    symmetryOptional = static_cast<Symmetry::Name>(
      temple::random::getSingle<unsigned>(0, 15, randomnessEngine())
    );

    std::geometric_distribution<unsigned> gd {0.2};
    assignmentOptional = gd(randomnessEngine());
  }

  // If a symmetry is specified, the bond number must match
  std::vector<molassembler::hashes::BondInformation> bonds;
  unsigned S;
  if(symmetryOptional) {
    S = Symmetry::size(*symmetryOptional);
  } else {
    S = temple::random::getSingle<unsigned>(1, 8, randomnessEngine());
  }

  for(unsigned i = 0; i < S; ++i) {
    bonds.emplace_back(genBondInformation());
  }
  std::sort(bonds.begin(), bonds.end());

  return {
    static_cast<Utils::ElementType>(
      temple::random::getSingle<unsigned>(1, 112, randomnessEngine())
    ),
    bonds,
    symmetryOptional,
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

    if(
      resultsMap.count(result) > 0
      && arguments != resultsMap.at(result)
    ) {
      BOOST_REQUIRE_MESSAGE(
        false,
        "Found overlapping result for different arguments to hashAtomEnvironment!\n"
      );
    }

    resultsMap.emplace(result, std::move(arguments));
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

    Molecule a = IO::read(currentFilePath.string());

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

    const Molecule m = IO::read(currentFilePath.string());
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
      std::string bPathString = base.string() + "_isomorphic.dot";
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

// Atom stereocenter assignments are part of strict equivalency
BOOST_AUTO_TEST_CASE(MoleculeBasicRSInequivalency) {
  // Build an asymmetric tetrahedral carbon
  Molecule a {Utils::ElementType::C, Utils::ElementType::H, BondType::Single};
  a.addAtom(Utils::ElementType::F, 0, BondType::Single);
  a.addAtom(Utils::ElementType::Cl, 0, BondType::Single);
  a.addAtom(Utils::ElementType::Br, 0, BondType::Single);
  a.setGeometryAtAtom(0, Symmetry::Name::Tetrahedral);

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

  // Set the geometries
  a.setGeometryAtAtom(0, Symmetry::Name::TrigonalPlanar);
  a.setGeometryAtAtom(1, Symmetry::Name::TrigonalPlanar);

  // Progression must recognize the new stereopermutator
  auto stereopermutatorOption = a.stereopermutators().option(
    BondIndex {0, 1}
  );
  BOOST_CHECK(
    stereopermutatorOption
    && stereopermutatorOption->numStereopermutations() == 2
  );

  std::cout << a << "\n";

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

  auto pseudocenter = IO::read(
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
   * to original abstract site case after forcing cut-tetrahedral symmetries
   * for each nitrogen.
   */
  auto complex = IO::read("various/propagation-test-case-1.json");

  for(AtomIndex i = 0; i < complex.graph().N(); ++i) {
    if(
      complex.graph().elementType(i) == Utils::ElementType::N
      && complex.stereopermutators().option(i)
      && Symmetry::size(complex.stereopermutators().option(i)->getSymmetry()) == 3
    ) {
      complex.setGeometryAtAtom(i, Symmetry::Name::CutTetrahedral);
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

  BOOST_REQUIRE_NO_THROW(molSplat = IO::split("multiple_molecules/multi_interpret.mol"));
  BOOST_REQUIRE_NO_THROW(xyzSplat = IO::split("multiple_molecules/multi_interpret.xyz"));

  // Each file contains one ethane, four water, and a potassium

  Molecule water {Utils::ElementType::O};
  water.addAtom(Utils::ElementType::H, 0, BondType::Single);
  water.addAtom(Utils::ElementType::H, 0, BondType::Single);
  water.setGeometryAtAtom(0, Symmetry::Name::Bent);

  Molecule ethane {Utils::ElementType::C};
  ethane.addAtom(Utils::ElementType::H, 0, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, 0, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, 0, BondType::Single);
  AtomIndex otherCarbon = ethane.addAtom(Utils::ElementType::C, 0, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, otherCarbon, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, otherCarbon, BondType::Single);
  ethane.addAtom(Utils::ElementType::H, otherCarbon, BondType::Single);
  ethane.setGeometryAtAtom(0, Symmetry::Name::Tetrahedral);
  ethane.setGeometryAtAtom(otherCarbon, Symmetry::Name::Tetrahedral);

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

  if(auto suggestedSymmetryOption = testMol.inferSymmetry(1u, stereocenterOption->getRanking())) {
    BOOST_CHECK(suggestedSymmetryOption.value() == Symmetry::Name::CutTetrahedral);
    testMol.setGeometryAtAtom(1u, suggestedSymmetryOption.value());
  }

  testMol.addAtom(Utils::ElementType::H, 1u, BondType::Single);

  BOOST_CHECK(
    testMol.stereopermutators().option(1u)
    && testMol.stereopermutators().option(1u)->getSymmetry() == Symmetry::Name::Tetrahedral
  );
}

// Isomer predicates work as expected
BOOST_AUTO_TEST_CASE(IsomerPredicateTests) {
  Molecule a, b;
  BOOST_REQUIRE_NO_THROW(a = molassembler::IO::read("isomers/enantiomers/Citalopram-R.mol"));
  BOOST_REQUIRE_NO_THROW(b = molassembler::IO::read("isomers/enantiomers/Citalopram-S.mol"));

  constexpr auto bitmask = AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Symmetries;

  a.canonicalize(bitmask);
  b.canonicalize(bitmask);

  BOOST_CHECK_MESSAGE(
    molassembler::enantiomeric(a, b),
    "Citalopram-R and -S are falsely determined not to be enantiomers."
  );

  BOOST_REQUIRE_NO_THROW(a = molassembler::IO::read("isomers/enantiomers/Isoleucine-RS.mol"));
  BOOST_REQUIRE_NO_THROW(b = molassembler::IO::read("isomers/enantiomers/Isoleucine-SR.mol"));

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
