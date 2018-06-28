#define BOOST_TEST_MODULE MoleculeTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/regex.hpp"

#include <iostream>

#include "detail/StdlibTypeAlgorithms.h"
#include "IO.h"

#include "temple/Random.h"
#include "temple/Invoke.h"
#include "temple/Stringify.h"

using namespace molassembler;

BOOST_AUTO_TEST_CASE( read_mol ) {
  std::vector<std::string> files {
    "test_files/2,2-dimethybutane.mol",
    "test_files/asymCarbon.mol",
    "test_files/C8H12_asym.mol",
    "test_files/opt-T-shaped0.mol",
    "test_files/opt-tetrahedral0.mol",
    "test_files/opt-trigonal-pyramidal0.mol",
    "test_files/opt-bent0.mol"
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
  std::vector<molassembler::BondType>,
  boost::optional<Symmetry::Name>,
  boost::optional<unsigned>
>;

HashArgumentsType randomArguments() {
  auto bonds = rng.getN<unsigned>(
    0,
    7,
    rng.getSingle<unsigned>(1, 8)
  );
  std::sort(bonds.begin(), bonds.end());

  boost::optional<Symmetry::Name> symmetryOptional;
  boost::optional<unsigned> assignmentOptional;
  if(rng.getSingle<bool>()) {
    symmetryOptional = static_cast<Symmetry::Name>(
      rng.getSingle<unsigned>(0, 15)
    );

    std::geometric_distribution<unsigned> gd {0.2};
    assignmentOptional = gd(rng.engine);
  }

  return {
    static_cast<Delib::ElementType>(
      rng.getSingle<unsigned>(1, 112)
    ),
    temple::cast<molassembler::BondType>(bonds),
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
  std::unordered_map<long long unsigned, HashArgumentsType> resultsMap;
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
          << "A: " << temple::stringify(arguments) << "\n"
          << "B: " << temple::stringify(resultsMap.at(result)) << "\n"
          << "both yield " << result << "\n"
      );
    }

    resultsMap.emplace(result, std::move(arguments));
  }
}

BOOST_AUTO_TEST_CASE(isomorphismTests) {
  const std::string directoryPrefix = "test_files/isomorphisms/"s;
  const boost::regex isomorphismFileRegex {R"(.+_isomorphism.mol)"};
  const boost::regex removeRegex {R"(_isomorphism.mol)"};

  /* TODO collect the original molecules in the folder and test all of them
   * against one another, ensuring !=
   */

  std::vector<Molecule> originals;

  boost::filesystem::path filesPath(directoryPrefix);
  boost::filesystem::recursive_directory_iterator end;
  for(boost::filesystem::recursive_directory_iterator iter(filesPath); iter != end; iter++) {
    const boost::filesystem::path currentFilePath = *iter;

    boost::smatch what;

    if(!boost::regex_match(iter -> path().filename().string(), what, isomorphismFileRegex)) {
      continue;
    }

    /* TODO Find some other way to compose the original string, i.e. regex replace?
     * replace the matching part with "" -> original path
     */

    auto originalFilePath = currentFilePath.parent_path() / (
      boost::regex_replace(iter -> path().filename().string(), removeRegex, "") + ".mol"
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

BOOST_AUTO_TEST_CASE(propagateGraphChangeTests) {
  const std::string directoryPrefix = "test_files/ranking_tree_molecules/"s;

  auto pseudocenter = IO::read(
    directoryPrefix + "(2R,3r,4S)-pentane-2,3,4-trithiol.mol"
  );

  AtomIndexType central = 0;
  std::array<AtomIndexType, 2> outer {{1, 6}};

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

  BOOST_CHECK(!pseudocenter.getStereocenterList().isStereogenic(central));

  // Make 6 from S to R -> stereocenter should reappear
  pseudocenter.assignStereocenter(outer.back(), 1);

  BOOST_CHECK(pseudocenter.getStereocenterList().isStereogenic(central));
}

BOOST_AUTO_TEST_CASE(moleculeSplitRecognition) {
  auto molSplat = IO::split("test_files/multiple_molecules/ethane_four_water.mol");

  auto xyzHandler = IO::FileHandlers::XYZHandler {};
  auto xyzData = xyzHandler.read("test_files/multiple_molecules/ethane_four_water.xyz");

  auto xyzSplat = IO::split("test_files/multiple_molecules/ethane_four_water.xyz");

  BOOST_CHECK(molSplat.size() == 5 && xyzSplat.size() == 5);
}
