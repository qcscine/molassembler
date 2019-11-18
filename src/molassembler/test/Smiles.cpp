/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "boost/test/unit_test.hpp"

#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/IO/SmilesParser.h"

#include "temple/Functional.h"
#include "temple/Stringify.h"

/* TODO
 * - allene "NC(Br)=[C@]=C(O)C"
     -> different smiles: {"F/C=C=C=C/F", R"(F/C=C=C=C\F)"}, // allene trans, cis
 * - Tests for
 *   - correct inferred shapes in simple cases
 *   - correct bond type and shape inferral in aromatic cycles and heterocycles
 *   - square planar centers
 *   - trigonal bipyramidal centers
 *   - octahedral centers
 */

using namespace Scine;
using namespace molassembler;

Molecule expectSingle(std::vector<Molecule>&& a) {
  if(a.size() == 1) {
    return a.front();
  }

  throw std::runtime_error("Expected single molecule result");
}

BOOST_AUTO_TEST_CASE(SmilesHydrogenFilling) {
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
    {"C#N", 3} // HCN
  };

  // TODO add more tests for P, S, N (these are the odd ones out): e.g. what is the right smiles for H2SO4?

  for(const auto& pair : pairs) {
    Molecule result;
    BOOST_REQUIRE_NO_THROW(result = expectSingle(IO::parseSmiles(pair.first)));
    BOOST_CHECK_MESSAGE(
      result.graph().N() ==  pair.second,
      "Expected " << pair.second << " atoms for '" << pair.first << "', got "
      << result.graph().N() << " instead."
    );
  }
}

BOOST_AUTO_TEST_CASE(SmilesClosesRingCycles) {
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
    BOOST_REQUIRE_NO_THROW(result = expectSingle(IO::parseSmiles(pair.first)));
    BOOST_CHECK_EQUAL(result.graph().B(), pair.second);
  }
}

BOOST_AUTO_TEST_CASE(AcceptValidSmiles) {
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
    BOOST_REQUIRE_NO_THROW(results = IO::parseSmiles(str));
    BOOST_CHECK(results.size() == 1);
  }
}

BOOST_AUTO_TEST_CASE(RejectInvalidSmiles) {
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
    BOOST_CHECK_THROW(IO::parseSmiles(str), std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(IdenticalSmiles) {
  const std::vector<std::pair<std::string, std::string>> pairs {
    {"C", "[CH4]"},
    {"[H][CH2][H]", "[H]C([H])([H])[H]"},
    {"C-C", "CC"},
    {"C-C-O", "CCO"},
    {"C-C=C-C", "CC=CC"},
    {"c1ccccc1", "C1=CC=CC=C1"},
    {"c1ccc2CCCc2c1", "C1=CC=CC(CCC2)=C12"},
    {"c1occc1", "C1OC=CC=1"},
    {"c1ccc1", "C1=CC=C1"},
    {"C1.C1", "CC"}, // Odd use of dot, but legal
    {"C1.C12.C2", "CCC"},
    {"c1c2c3c4cc1.Br2.Cl3.Cl4", "C1=CC(=C(C(=C1)Br)Cl)Cl"},
    {"N[C@](Br)(O)C", "Br[C@](O)(N)C"},
    {"O[C@](Br)(C)N", "Br[C@](C)(O)N"},
    {"C[C@](Br)(N)O", "Br[C@](N)(C)O"},
    {"C[C@@](Br)(O)N", "Br[C@@](N)(O)C"},
    {"[C@@](C)(Br)(O)N", "[C@@](Br)(N)(O)C"},
    {"FC1C[C@](Br)(Cl)CCC1", "C@]1(Br)(Cl)CCCC(F)C1"},
    {"F/C=C/F", R"y(F\C=C\F)y"},
    {"F/C=C/F", R"y(C(\F)=C/F)y"},
    {R"y(F\C=C/F)y", R"y(C(/F)=C/F)y"},
  };

  for(const auto& pair : pairs) {
    Molecule a, b;
    BOOST_REQUIRE_NO_THROW(a = expectSingle(IO::parseSmiles(pair.first)));
    BOOST_REQUIRE_NO_THROW(b = expectSingle(IO::parseSmiles(pair.second)));
    BOOST_CHECK_MESSAGE(
      a == b,
      "Smiles pair " << pair.first << ", " << pair.second << " did not compare equal as expected"
    );
  }
}

BOOST_AUTO_TEST_CASE(DifferentSmiles) {
  // Pairs of smiles that give different molecules
  const std::vector<std::pair<std::string, std::string>> pairs {
    {"F/C=C/F", R"y(C(/F)=C/F)y"}, // trans, cis
    {"F/C=C/F", R"y(F\C=C/F)y"}, // trans, cis
    {R"y(C(\F)=C/F)y", R"y(F\C=C/F)y"}, // trans, cis
    {R"y(C(\F)=C/F)y", R"y(C(/F)=C/F)y"}, // trans, cis
  };

  for(const auto& pair : pairs) {
    Molecule a, b;
    BOOST_REQUIRE_NO_THROW(a = expectSingle(IO::parseSmiles(pair.first)));
    BOOST_REQUIRE_NO_THROW(b = expectSingle(IO::parseSmiles(pair.second)));
    BOOST_CHECK_MESSAGE(
      a != b,
      "Smiles pair " << pair.first << ", " << pair.second << " did not compare different as expected"
    );
    std::cout << "A: " << a << "\nB:" << b << "\n\n";
  }
}

BOOST_AUTO_TEST_CASE(SmilesWithMultipleMolecules) {
  const std::vector<std::pair<std::string, std::vector<unsigned>>> pairs {
    {"[Na+].[Cl-]", {{1, 1}}},
    {"Oc1ccccc1.NCCO", {{11, 13}}}, // Regular spec of multiple molecules
    {"c1cc(O.NCCO)ccc1", {{11, 13}}}, // Irregular, but valid
    {"Oc1cc(.NCCO)ccc1", {{11, 13}}}, // even more irregular, but valid
    {"[NH4+].[NH4+].[O-]S(=O)(=O)[S-]", {{5, 5, 5}}},
    {"C=1CCCCC1.C=1CCCCC1", {{16, 16}}} // Reuse of ring closure numbers
  };

  // parse, count sizes, order and lex. compare
  for(const auto& pair : pairs) {
    std::vector<Molecule> results;
    BOOST_REQUIRE_NO_THROW(results = IO::parseSmiles(pair.first));
    BOOST_REQUIRE(results.size() > 1);

    BOOST_CHECK_EQUAL(results.size(), pair.second.size());

    auto sizes = temple::sort(
      temple::map(
        results,
        [](const Molecule& m) -> unsigned { return m.graph().N(); }
      )
    );

    BOOST_CHECK_MESSAGE(
      sizes == pair.second,
      "Expected sizes: " << temple::stringify(pair.second) << ", got " << temple::stringify(sizes) << " instead for smiles '" << pair.first << "'"
    );
  }
}
