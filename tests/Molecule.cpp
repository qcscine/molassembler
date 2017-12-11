#define BOOST_TEST_MODULE MoleculeTests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include "IO.h"
#include "StdlibTypeAlgorithms.h"

BOOST_AUTO_TEST_CASE( read_mol ) {
  using namespace MoleculeManip;

  std::vector<std::string> files {
    "../tests/mol_files/2,2-dimethybutane.mol",
    "../tests/mol_files/asymCarbon.mol",
    "../tests/mol_files/C8H12_asym.mol",
    "../tests/mol_files/opt-T-shaped0.mol",
    "../tests/mol_files/opt-tetrahedral0.mol",
    "../tests/mol_files/opt-trigonal-pyramidal0.mol",
    "../tests/mol_files/opt-bent0.mol"
  };

  // instantiate reader
  IO::MOLFileHandler molHandler;

  for(const auto& filename : files) {
    Molecule mol = molHandler.readSingle(filename);
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
  // Molecule should be trivially usable on the stack
  using namespace MoleculeManip;

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
    TemplateMagic::all_of(
      TemplateMagic::mapAllPairs(
        allJustHydrogen,
        std::equal_to<Molecule>()
      )
    ),
    "Basic construction does not generate equal molecules"
  );
}

const std::string directoryPrefix = "../tests/mol_files/ranking_tree_molecules/"s;

BOOST_AUTO_TEST_CASE(propagateGraphChangeTests) {
  using namespace MoleculeManip;
  IO::MOLFileHandler molReader;

  auto pseudocenter = molReader.readSingle(
    directoryPrefix + "(2R,3r,4S)-pentane-2,3,4-trithiol.mol"
  );

  AtomIndexType central = 0;
  std::array<AtomIndexType, 2> outer {1, 6};

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
  pseudocenter.assignStereocenterAtAtom(outer.front(), 0);

  BOOST_CHECK(!pseudocenter.getStereocenterList().involving(central));

  // Make 6 from S to R -> stereocenter should reappear
  pseudocenter.assignStereocenterAtAtom(outer.back(), 1);

  BOOST_CHECK(pseudocenter.getStereocenterList().involving(central));
}
