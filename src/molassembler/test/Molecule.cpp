#define BOOST_TEST_MODULE MoleculeTests
#define BOOST_FILESYSTEM_NO_DEPRECATED

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/test/unit_test.hpp"

#include "chemical_symmetries/Symmetries.h"

#include "temple/Containers.h"
#include "temple/Random.h"
#include "temple/Invoke.h"
#include "temple/Stringify.h"

#include "molassembler/Detail/StdlibTypeAlgorithms.h"
#include "molassembler/IO/FileHandlers.h"
#include "molassembler/Molecule.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/Options.h"
#include "molassembler/StereocenterList.h"

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
      temple::mapAllPairs(
        allJustHydrogen,
        std::equal_to<>()
      )
    ),
    "Basic construction does not generate equal molecules"
  );
}

using HashArgumentsType = std::tuple<
  Delib::ElementType,
  std::vector<molassembler::hashes::BondInformation>,
  boost::optional<Symmetry::Name>,
  boost::optional<unsigned>
>;

HashArgumentsType randomArguments() {
  auto genBondInformation = []() -> molassembler::hashes::BondInformation {
    BondType bty = static_cast<BondType>(prng.getSingle<unsigned>(0, 6));

    bool bondStereocenterPresent = prng.getSingle<bool>();
    boost::optional<unsigned> bondStereocenterAssignment = boost::none;

    if(bondStereocenterPresent) {
      bondStereocenterAssignment = prng.getSingle<unsigned>(0, 1);
    }

    return {
      bty,
      bondStereocenterPresent,
      bondStereocenterAssignment
    };
  };

  boost::optional<Symmetry::Name> symmetryOptional;
  boost::optional<unsigned> assignmentOptional;
  if(prng.getSingle<bool>()) {
    symmetryOptional = static_cast<Symmetry::Name>(
      prng.getSingle<unsigned>(0, 15)
    );

    std::geometric_distribution<unsigned> gd {0.2};
    assignmentOptional = gd(prng.engine);
  }

  // If a symmetry is specified, the bond number must match
  std::vector<molassembler::hashes::BondInformation> bonds;
  unsigned S;
  if(symmetryOptional) {
    S = Symmetry::size(*symmetryOptional);
  } else {
    S = prng.getSingle<unsigned>(1, 8);
  }

  for(unsigned i = 0; i < S; ++i) {
    bonds.emplace_back(genBondInformation());
  }
  std::sort(bonds.begin(), bonds.end());

  return {
    static_cast<Delib::ElementType>(
      prng.getSingle<unsigned>(1, 112)
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

  /* TODO collect the original molecules in the folder and test all of them
   * against one another, ensuring !=
   */

  std::vector<Molecule> originals;

  for(
    const boost::filesystem::path& currentFilePath :
    boost::filesystem::recursive_directory_iterator(directoryBase)
  ) {
    boost::smatch what;

    if(!boost::regex_match(currentFilePath.filename().string(), what, isomorphismFileRegex)) {
      continue;
    }

    /* TODO Find some other way to compose the original string, i.e. regex replace?
     * replace the matching part with "" -> original path
     */

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
      temple::mapAllPairs(
        originals,
        std::not_equal_to<> {}
      )
    ),
    "Some originals in the isomorphism test folder match one another!"
  );
}

BOOST_AUTO_TEST_CASE(basicRSInequivalencyTests) {
  // Build an asymmetric tetrahedral carbon
  Molecule a {Delib::ElementType::C, Delib::ElementType::H, BondType::Single};
  a.addAtom(Delib::ElementType::F, 0, BondType::Single);
  a.addAtom(Delib::ElementType::Cl, 0, BondType::Single);
  a.addAtom(Delib::ElementType::Br, 0, BondType::Single);
  a.setGeometryAtAtom(0, Symmetry::Name::Tetrahedral);

  // Make sure it's recognized as asymmetric
  auto centralStereocenterOption = a.stereocenters().option(0);
  BOOST_CHECK(
    centralStereocenterOption
    && centralStereocenterOption->numStereopermutations() == 2
  );

  // Assign it and create its opposite stereopermutation in another Molecule
  a.assignStereocenter(0, 0);
  Molecule b = a;
  b.assignStereocenter(0, 1);

  // These must compare unequal
  BOOST_CHECK(a != b);
}

BOOST_AUTO_TEST_CASE(basicEZInequivalencyTests) {
  Molecule a {Delib::ElementType::C, Delib::ElementType::C, BondType::Double};
  a.addAtom(Delib::ElementType::H, 0, BondType::Single);
  a.addAtom(Delib::ElementType::F, 0, BondType::Single);
  a.addAtom(Delib::ElementType::H, 1, BondType::Single);
  a.addAtom(Delib::ElementType::F, 1, BondType::Single);

  // Set the geometries
  a.setGeometryAtAtom(0, Symmetry::Name::TrigonalPlanar);
  a.setGeometryAtAtom(1, Symmetry::Name::TrigonalPlanar);

  // Progression must recognize the new stereocenter
  auto stereocenterOption = a.stereocenters().option(
    BondIndex {0, 1}
  );
  BOOST_CHECK(
    stereocenterOption
    && stereocenterOption->numStereopermutations() == 2
  );

  std::cout << a << "\n";

  a.assignStereocenter(
    BondIndex {0, 1},
    0
  );

  Molecule b = a;
  b.assignStereocenter(
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
  auto stereocenterOption = molecule.stereocenters().option(i);

  if(!stereocenterOption) {
    return false;
  }

  if(stereocenterOption->numStereopermutations() <= 1) {
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

  /* If the outer stereocenters have the same assignment, the central
   * stereocenter shouldn't exist. If the outer stereocenters have a different
   * assignment, the central stereocenter has to exist.
   *
   * Initially, we have
   * 0 -> 1 r
   * 1 -> 1 R
   * 6 -> 0 S
   */
  // Make 1 from R to S -> stereocenter should disappear
  pseudocenter.assignStereocenter(outer.front(), 0);

  BOOST_CHECK(!isStereogenic(pseudocenter, central));

  // Make 6 from S to R -> stereocenter should reappear
  pseudocenter.assignStereocenter(outer.back(), 1);

  BOOST_CHECK(isStereogenic(pseudocenter, central));
}

BOOST_AUTO_TEST_CASE(moleculeSplitRecognition) {
  auto molSplat = IO::split("multiple_molecules/ethane_four_water.mol");

  auto xyzHandler = IO::FileHandlers::XYZHandler {};
  auto xyzData = xyzHandler.read("multiple_molecules/ethane_four_water.xyz");

  auto xyzSplat = IO::split("multiple_molecules/ethane_four_water.xyz");

  BOOST_CHECK(molSplat.size() == 5 && xyzSplat.size() == 5);
}
