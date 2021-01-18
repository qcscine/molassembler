/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "boost/test/unit_test.hpp"

#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"
#include "Molassembler/IO/SmilesParser.h"
#include "Molassembler/IO/SmilesEmitter.h"
#include "Molassembler/IO.h"

#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Stringify.h"

#include "Fixtures.h"

#include <iostream>

/* TODO
 * - allene "NC(Br)=[C@]=C(O)C"
     -> different smiles: {"F/C=C=C=C/F", R"(F/C=C=C=C\F)"}, // allene trans, cis
 * - Tests for
 *   - correct inferred shapes in simple cases
 *   - correct bond type and shape inferral in aromatic cycles and heterocycles
 */

using namespace Scine::Molassembler;

Molecule expectSingle(std::vector<Molecule>&& a) {
  if(a.size() == 1) {
    return a.front();
  }

  throw std::runtime_error("Expected single molecule result");
}

BOOST_AUTO_TEST_CASE(SmilesHydrogenFilling, *boost::unit_test::label("Molassembler")) {
  // Pairs of smiles strings and atom counts
  std::vector<std::pair<std::string, unsigned>> pairs {
    {"C", 5},
    {"[CH4]", 5},
    {"N", 4},
    {"[NH3]", 4},
    {"Cl", 2},
    {"[ClH]", 2},
    {"[H]Cl", 2},
    {"[H][Cl]", 2},
    {"BrBr", 2},
    {"FCl", 2},
    {"Br", 2},
    {"I", 2},
    {"F", 2},
    {"[F]", 1}, // Atom brackets aren't valence filled
    {"CC", 8}, // C2H6
    {"C=C", 6}, // C2H4
    {"C#C", 4}, // C2H2
    {"C=O", 4}, // H2CO
    {"C#N", 3}, // HCN
    {"N(C)(C)(C)C", 18}, // HN(Me)4 (tricky one for valence either 3 or 5)
  };

  // TODO add more tests for P, S, N (these are the odd ones out): e.g. what is the right smiles for H2SO4?

  for(const auto& pair : pairs) {
    Molecule result;
    BOOST_REQUIRE_NO_THROW(result = expectSingle(IO::Experimental::parseSmiles(pair.first)));
    BOOST_CHECK_MESSAGE(
      result.graph().V() ==  pair.second,
      "Expected " << pair.second << " atoms for '" << pair.first << "', got "
      << result.graph().V() << " instead."
    );
  }
}

BOOST_AUTO_TEST_CASE(SmilesClosesRingCycles, *boost::unit_test::label("Molassembler")) {
  // Pairs of Smiles strings and expected numbers of edges
  std::vector<std::pair<std::string, unsigned>> pairs {
    {"C1=CC=CC=C1", 12}, // benzene
    {"C1CCCCC=1", 16}, // three variations on cyclohexene, all valid
    {"C=1CCCCC1", 16},
    {"C=1CCCCC=1", 16}
  };

  // TODO more ring-closing tests

  for(const auto& pair : pairs) {
    Molecule result;
    BOOST_REQUIRE_NO_THROW(result = expectSingle(IO::Experimental::parseSmiles(pair.first)));
    BOOST_CHECK_EQUAL(result.graph().E(), pair.second);
  }
}

BOOST_AUTO_TEST_CASE(AcceptValidSmiles, *boost::unit_test::label("Molassembler")) {
  const std::vector<std::string> validSmiles {
    "[HH0]",
    "[H][H]",
    "C",
    "CC",
    "CCO",
    "NCCCC",
    "CCCCN",
    "C=C",
    "C#N",
    "CC#CC",
    "CCC=O",
    "[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl",
    "CCC(CC)CO",
    "CC(C)C(=O)C(C)C",
    "OCC(CCC)C(C(C)C)CCC",
    "OS(=O)(=S)O",
    "C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C",
    "C1CCCCC1",
    "N1CC2CCCC2CC1",
    "C=1CCCCC=1",
    "C=1CCCCC1",
    "C1CCCCC=1",
    "C1CCCCC1C1CCCCC1",
    "C12(CCCCC1)CCCCC2",
    "F/C(CC)=C/F",
  };

  std::vector<Molecule> results;
  for(const auto& str : validSmiles) {
    BOOST_REQUIRE_NO_THROW(results = IO::Experimental::parseSmiles(str));
    BOOST_CHECK(results.size() == 1);
  }
}

BOOST_AUTO_TEST_CASE(RejectInvalidSmiles, *boost::unit_test::label("Molassembler")) {
  const std::vector<std::string> invalidSmiles {
    "[HH]", // hydrogens atoms cannot have hydrogen counts
    "[HH1]",
    "[HH3]",
    "C1CCC", // unmatched ring closure
    "C-1CCCCCC=1", // mismatched explicit ring closure bond types
    "C12CCCC12", // parallel ring closure bond (2x same atoms)
    "C12C2CCC1", // parallel bonds (existing implicit + ring closure)
    "C11", // self-loop ring closure
    "CC(CC", // unmatched branch open
    "CC)CC", // unmatched branch closure
    "F.FB(F)F.[NH2+251][C@@H](CP(c1ccccc1)c1ccccc1)C(C)(C)C", // odd trailing characters in atom bracket
    "C.1CCCCC.1", // dot before ring closing bond
    R"y(C/C(\F)=C/F)y", // both C and F are "down" at the left carbon
  };

  for(const auto& str : invalidSmiles) {
    BOOST_CHECK_THROW(IO::Experimental::parseSmiles(str), std::runtime_error);
  }
}

/* This needs to be executed in the low temperature regime for the trigonal
 * bipyramidal smiles (these are considered achiral in the high temperature
 * approximation)
 */
BOOST_FIXTURE_TEST_CASE(IdenticalSmiles, LowTemperatureFixture) {
  const std::vector<std::pair<std::string, std::string>> pairs {
    {"C", "[CH4]"}, // Implicit hydrogens
    {"[H][CH2][H]", "[H]C([H])([H])[H]"},
    {"C-C", "CC"}, // Implicit bonds
    {"C-C-O", "CCO"},
    {"C-C=C-C", "CC=CC"},
    {"C1.C1", "CC"}, // Odd use of dot, but legal
    {"C1.C12.C2", "CCC"},
    {"N[C@](Br)(O)C", "Br[C@](O)(N)C"}, // Tetrahedral stereo
    {"O[C@](Br)(C)N", "Br[C@](C)(O)N"},
    {"C[C@](Br)(N)O", "Br[C@](N)(C)O"},
    {"C[C@@](Br)(O)N", "Br[C@@](N)(O)C"},
    {"[C@@](C)(Br)(O)N", "[C@@](Br)(N)(O)C"},
    {"N[C@H](O)C", "[H][C@](N)(O)C"},
    {"F/C=C/F", R"y(F\C=C\F)y"}, // Double bond stereo
    {"F/C=C/F", R"y(C(\F)=C/F)y"},
    {R"y(F\C=C/F)y", R"y(C(/F)=C/F)y"},
    {"S[As@TB1](F)(Cl)(Br)N", "S[As@TB2](Br)(Cl)(F)N"}, // Trig bipy stereo
    {"S[As@TB5](F)(N)(Cl)Br", "F[As@TB10](S)(Cl)(N)Br"},
    {"F[As@TB15](Cl)(S)(Br)N", "Br[As@TB20](Cl)(S)(F)N"},
    {"C[Co@](F)(Cl)(Br)(I)S", "F[Co@@](S)(I)(C)(Cl)Br"}, // Octahedral stereo
    {"S[Co@OH5](F)(I)(Cl)(C)Br", "Br[Co@OH9](C)(S)(Cl)(F)I"},
    {"Br[Co@OH12](Cl)(I)(F)(S)C", "Cl[Co@OH15](C)(Br)(F)(I)S"},
    {"Cl[Co@OH19](C)(I)(F)(S)Br", "I[Co@OH27](Cl)(Br)(F)(S)C"},
  };

  // TODO known failures
  //   {"FC1C[C@](Br)(Cl)CCC1", "[C@]1(Br)(Cl)CCCC(F)C1"}, // Peculiarity with ring closing bond stereo order
  //   {R"y(F/C=C/1.Br1)y", R"y(F/C=C1.Br\1)y"}, // Funky use of dot
  //   {"[C@@](O)(Cl)(F)Br", "[C@@](O)(Cl)(F)1.Br1"}, // Funky use of dot
  for(const auto& pair : pairs) {
    Molecule a, b;
    BOOST_REQUIRE_NO_THROW(a = expectSingle(IO::Experimental::parseSmiles(pair.first)));
    BOOST_REQUIRE_NO_THROW(b = expectSingle(IO::Experimental::parseSmiles(pair.second)));
    BOOST_CHECK_MESSAGE(
      a == b,
      "Smiles pair " << pair.first << ", " << pair.second << " did not compare equal as expected"
    );

    if(a != b) {
      std::cout << "A: " << a << "\nB:" << b << "\n\n";
    }
  }
}

BOOST_AUTO_TEST_CASE(SimilarSmiles, *boost::unit_test::label("Molassembler")) {
  const std::vector<std::pair<std::string, std::string>> pairs {
    {"c1ccccc1", "C1=CC=CC=C1"}, // Benzene
    {"c1ccc2CCCc2c1", "C1=CC=CC(CCC2)=C12"}, // Indane
    {"c1occc1", "C1OC=CC=1"}, // Furan
    {"c1cncnc1", "C1=CN=CN=C1"}, // Pyrimidine
    {"c1cnccc1", "C1=CN=CC=C1"}, // Pyridine
    {"c1ccpcc1", "C1=CC=PC=C1"}, // Phosphorine
    {"n1c3c(cc2c1cccc2)cccc3", "C1=CC=C2C(=C1)C=C3C=CC=CC3=N2"}, // Acridine
    {"o2cc1ccccc1c2", "C1=CC2=COC=C2C=C1"}, // Isobenzofuran
    {"n2oc1ccccc1c2", "C1=CC=C2C(=C1)C=NO2"}, // Benzisoxazole
    {"c1cocn1", "C1=COC=N1"}, // Oxazole
    {"c1ccc1", "C1=CC=C1"}, // Cyclobutadiene (Not aromatic but valid SMILES)
    {"c1c2c3c4cc1.Br2.Cl3.Cl4", "C1=CC(=C(C(=C1)Br)Cl)Cl"}, // Aromatics + dot
  };

  const AtomEnvironmentComponents strictness = (
    AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::Shapes
  );

  for(const auto& pair : pairs) {
    Molecule a, b;
    BOOST_REQUIRE_NO_THROW(a = expectSingle(IO::Experimental::parseSmiles(pair.first)));
    BOOST_REQUIRE_NO_THROW(b = expectSingle(IO::Experimental::parseSmiles(pair.second)));

    const bool closeEnough = static_cast<bool>(a.modularIsomorphism(b, strictness));
    BOOST_CHECK_MESSAGE(
      closeEnough,
      "Smiles pair " << pair.first << ", " << pair.second << " did not have same connectivity and shapes as expected"
    );

    if(!closeEnough) {
      std::cout << "A: " << a << "\nB:" << b << "\n\n";
    }
  }
}

BOOST_AUTO_TEST_CASE(DifferentSmiles, *boost::unit_test::label("Molassembler")) {
  // Pairs of smiles that give different molecules
  const std::vector<std::pair<std::string, std::string>> pairs {
    {"F/C=C/F", R"y(C(/F)=C/F)y"}, // trans, cis
    {"F/C=C/F", R"y(F\C=C/F)y"}, // trans, cis
    {R"y(C(\F)=C/F)y", R"y(F\C=C/F)y"}, // trans, cis
    {R"y(C(\F)=C/F)y", R"y(C(/F)=C/F)y"}, // trans, cis
    {"N[C@](Br)(O)C", "N[C@@](Br)(O)C"}, // R, S
  };

  for(const auto& pair : pairs) {
    Molecule a;
    Molecule b;
    BOOST_REQUIRE_NO_THROW(a = expectSingle(IO::Experimental::parseSmiles(pair.first)));
    BOOST_REQUIRE_NO_THROW(b = expectSingle(IO::Experimental::parseSmiles(pair.second)));
    BOOST_CHECK_MESSAGE(
      a != b,
      "Smiles pair " << pair.first << ", " << pair.second << " did not compare different as expected"
    );
  }
}

BOOST_AUTO_TEST_CASE(SmilesWithMultipleMolecules, *boost::unit_test::label("Molassembler")) {
  const std::vector<std::pair<std::string, std::vector<unsigned>>> pairs {
    {"[Na+].[Cl-]", {{1, 1}}},
    {"[NH4+].[NH4+].[O-]S(=O)(=O)[S-]", {{5, 5, 5}}},
    {"C=1CCCCC1.C=1CCCCC1", {{16, 16}}}, // Reuse of ring closure numbers
    {"Oc1ccccc1.NCCO", {{11, 13}}}, // Regular spec of multiple molecules
    {"c1cc(O.NCCO)ccc1", {{11, 13}}}, // Irregular, but valid
    {"Oc1cc(.NCCO)ccc1", {{11, 13}}}, // even more irregular, but valid
  };

  // parse, count sizes, order and lex. compare
  for(const auto& pair : pairs) {
    std::vector<Molecule> results;
    BOOST_REQUIRE_NO_THROW(results = IO::Experimental::parseSmiles(pair.first));
    BOOST_REQUIRE(results.size() > 1);

    BOOST_CHECK_EQUAL(results.size(), pair.second.size());

    const auto sizes = Temple::sorted(
      Temple::map(
        results,
        [](const Molecule& m) -> unsigned { return m.graph().V(); }
      )
    );

    BOOST_CHECK_MESSAGE(
      sizes == pair.second,
      "Expected sizes: " << Temple::stringify(pair.second) << ", got " << Temple::stringify(sizes) << " instead for smiles '" << pair.first << "'"
    );
  }
}

BOOST_AUTO_TEST_CASE(EmitAliphatic, *boost::unit_test::label("Molassembler")) {
  const std::vector<std::string> cases {
    "[H][H]",
    "C",
    "CC",
    "CC(C)(C)C",
    "CN",
    "CNO",
    "C1CC1",
    "C1CCC1",
    "C1CCCC1"
  };

  for(const std::string& smiles : cases) {
    Molecule mol;
    Molecule mol2;
    std::string emitted;
    BOOST_REQUIRE_NO_THROW(mol = expectSingle(IO::Experimental::parseSmiles(smiles)));
    BOOST_REQUIRE_NO_THROW(emitted = IO::Experimental::emitSmiles(mol));
    std::cout << smiles << " -> " << emitted << "\n";
    BOOST_REQUIRE_NO_THROW(mol2 = expectSingle(IO::Experimental::parseSmiles(emitted)));
    BOOST_CHECK(mol == mol2);
  }
}
