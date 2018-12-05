// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#define BOOST_TEST_MODULE MoleculeTests
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/test/unit_test.hpp"

#include "chemical_symmetries/Symmetries.h"

#include "temple/Adaptors/AllPairs.h"
#include "temple/Functional.h"
#include "temple/Random.h"
#include "temple/Invoke.h"
#include "temple/Stringify.h"

#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/IO/FileHandlers.h"
#include "molassembler/Isomers.h"
#include "molassembler/Molecule.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/Options.h"
#include "molassembler/StereopermutatorList.h"

#include <iostream>

using namespace molassembler;

BOOST_AUTO_TEST_CASE( read_mol ) {
  std::vector<std::string> files {
    "2,2-dimethybutane.mol",
    "asymCarbon.mol",
    "C8H12_asym.mol",
    "opt-T-shaped0.mol",
    "opt-tetrahedral0.mol",
    "opt-trigonal-pyramidal0.mol",
    "opt-bent0.mol"
  };

  for(const auto& filename : files) {
    Molecule mol = IO::read(filename);
    // Invoke ostream operator
    std::cout << mol << std::endl;

    // Make dot files for every file
    auto slashSplat = StdlibTypeAlgorithms::split(filename, '/');
    auto dotSplat = StdlibTypeAlgorithms::split(slashSplat.back(), '.');

    std::ofstream dotFile (dotSplat.front() + ".dot");
    dotFile << mol.dumpGraphviz();
    dotFile.close();
  }
}

BOOST_AUTO_TEST_CASE(ruleOfFiveTrivial) {

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
  Scine::Utils::ElementType,
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
    static_cast<Scine::Utils::ElementType>(
      temple::random::getSingle<unsigned>(1, 112, randomnessEngine())
    ),
    bonds,
    symmetryOptional,
    assignmentOptional
  };
}

BOOST_AUTO_TEST_CASE(environmentHashingTests) {
  auto bitmaskTuple = std::make_tuple(
    temple::make_bitmask(molassembler::AtomEnvironmentComponents::ElementTypes)
      | molassembler::AtomEnvironmentComponents::BondOrders
      | molassembler::AtomEnvironmentComponents::Symmetries
      | molassembler::AtomEnvironmentComponents::Stereopermutations
  );

  // Try to guess a disjoint combination that has the same value
  std::unordered_map<molassembler::hashes::WideHashType, HashArgumentsType> resultsMap;
  for(unsigned N = 0; N < 1e6; ++N) {
    auto arguments = randomArguments();

    auto result = temple::detail::invokeHelper(
      hashes::atomEnvironment,
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

BOOST_AUTO_TEST_CASE(isomorphismTests) {
  using namespace std::string_literals;

  boost::filesystem::path directoryBase("isomorphisms");

  const boost::regex isomorphismFileRegex {R"(.+_isomorphism.mol)"};
  const boost::regex removeRegex {R"(_isomorphism.mol)"};

  std::vector<Molecule> originals;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    boost::smatch what;

    if(!boost::regex_match(currentFilePath.filename().string(), what, isomorphismFileRegex)) {
      continue;
    }

    auto originalFilePath = currentFilePath.parent_path() / (
      boost::regex_replace(currentFilePath.filename().string(), removeRegex, "") + ".mol"
    );

    if(!boost::filesystem::exists(originalFilePath)) {
      std::cout << "There is no matching file " << originalFilePath
        << " to " << currentFilePath << std::endl;
      continue;
    }

    Molecule a = IO::read(originalFilePath.string());
    Molecule b = IO::read(currentFilePath.string());

    BOOST_CHECK_MESSAGE(
      a == b,
      "Molecule isomorphism fails for " << what[0].str() << "!"
    );

    originals.push_back(std::move(a));
  }

  BOOST_CHECK_MESSAGE(
    temple::all_of(
      temple::adaptors::allPairs(originals),
      std::not_equal_to<> {}
    ),
    "Some originals in the isomorphism test folder match one another!"
  );
}

BOOST_AUTO_TEST_CASE(basicRSInequivalencyTests) {
  // Build an asymmetric tetrahedral carbon
  Molecule a {Scine::Utils::ElementType::C, Scine::Utils::ElementType::H, BondType::Single};
  a.addAtom(Scine::Utils::ElementType::F, 0, BondType::Single);
  a.addAtom(Scine::Utils::ElementType::Cl, 0, BondType::Single);
  a.addAtom(Scine::Utils::ElementType::Br, 0, BondType::Single);
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

BOOST_AUTO_TEST_CASE(basicEZInequivalencyTests) {
  Molecule a {Scine::Utils::ElementType::C, Scine::Utils::ElementType::C, BondType::Double};
  a.addAtom(Scine::Utils::ElementType::H, 0, BondType::Single);
  a.addAtom(Scine::Utils::ElementType::F, 0, BondType::Single);
  a.addAtom(Scine::Utils::ElementType::H, 1, BondType::Single);
  a.addAtom(Scine::Utils::ElementType::F, 1, BondType::Single);

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

bool isStereogenic(
  const Molecule& molecule,
  AtomIndex i
) {
  auto stereopermutatorOption = molecule.stereopermutators().option(i);

  if(!stereopermutatorOption) {
    return false;
  }

  if(stereopermutatorOption->numStereopermutations() <= 1) {
    return false;
  }

  return true;
}

BOOST_AUTO_TEST_CASE(propagateGraphChangeTests) {
  boost::filesystem::path filePath("ranking_tree_molecules");
  filePath /= "(2R,3r,4S)-pentane-2,3,4-trithiol.mol";

  auto pseudocenter = IO::read(
    filePath.string()
  );

  AtomIndex central = 0;
  std::array<AtomIndex, 2> outer {{1, 6}};

  /* If the outer stereopermutators have the same assignment, the central
   * stereopermutator shouldn't exist. If the outer stereopermutators have a different
   * assignment, the central stereopermutator has to exist.
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
}

BOOST_AUTO_TEST_CASE(moleculeSplitRecognition) {
  auto molSplat = IO::split("multiple_molecules/ethane_four_water.mol");

  auto xyzHandler = IO::FileHandlers::XYZHandler {};
  auto xyzData = xyzHandler.read("multiple_molecules/ethane_four_water.xyz");

  auto xyzSplat = IO::split("multiple_molecules/ethane_four_water.xyz");

  BOOST_CHECK(molSplat.size() == 5 && xyzSplat.size() == 5);
}

BOOST_AUTO_TEST_CASE(moleculeGeometryChoices) {
  molassembler::Molecule testMol(Scine::Utils::ElementType::Ru, Scine::Utils::ElementType::N, BondType::Single);
  testMol.addAtom(Scine::Utils::ElementType::H, 1u, BondType::Single);
  testMol.addAtom(Scine::Utils::ElementType::H, 1u, BondType::Single);

  BOOST_CHECK(
    testMol.stereopermutators().option(1u)
    && testMol.stereopermutators().option(1u)->getSymmetry() == Symmetry::Name::CutTetrahedral
  );

  testMol.addAtom(Scine::Utils::ElementType::H, 1u, BondType::Single);

  BOOST_CHECK(
    testMol.stereopermutators().option(1u)
    && testMol.stereopermutators().option(1u)->getSymmetry() == Symmetry::Name::Tetrahedral
  );
}

BOOST_AUTO_TEST_CASE(IsomerPredicateTests) {
  BOOST_CHECK_MESSAGE(
    molassembler::enantiomeric(
      molassembler::IO::read("isomers/enantiomers/Citalopram-R.mol"),
      molassembler::IO::read("isomers/enantiomers/Citalopram-S.mol"),
      molassembler::SameIndexingTag {}
    ),
    "Citalopram-R and -S are falsely determined not to be enantiomers."
  );

  BOOST_CHECK_MESSAGE(
    molassembler::enantiomeric(
      molassembler::IO::read("isomers/enantiomers/Isoleucine-RS.mol"),
      molassembler::IO::read("isomers/enantiomers/Isoleucine-SR.mol"),
      molassembler::SameIndexingTag {}
    ),
    "Citalopram-R and -S are falsely determined not to be enantiomers."
  );
}
