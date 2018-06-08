#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>
#include <numeric>

#include "stereopermutation/GenerateUniques.h"
#include "stereopermutation/LogicalOperatorTests.h"

#include "temple/Random.h"
#include "temple/Containers.h"

using namespace stereopermutation;

/* TODO
 * - add more tests from different geometries
 */

template<typename T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& vector) {
  for(unsigned i = 0; i < vector.size(); i++) {
    os << vector[i];
    if(i != vector.size() - 1) {
      os << ", ";
    }
  }

  return os;
}

// create instances of all symmetries with monodentate ligands
BOOST_AUTO_TEST_CASE( assignmentInstantiation ) {
  for(const auto& symmetryName: Symmetry::allNames) {
    Stereopermutation testStereopermutation(
      symmetryName,
      std::vector<char>(
        Symmetry::size(symmetryName),
        'A'
      )
    );
  }
}

BOOST_AUTO_TEST_CASE( assignment_basics ) {
  // Constructors
  Stereopermutation instanceWithBondedLigands(
    Symmetry::Name::Octahedral,
    {'A', 'B', 'C', 'D', 'E', 'F'},
    {
      std::make_pair(0,1),
      std::make_pair(2,3),
      std::make_pair(4,5)
    }
  );

  { // columnSwap
    auto instanceCopy = instanceWithBondedLigands;
    instanceCopy.columnSwap(0, 1);
    BOOST_CHECK(
      instanceCopy.toString()
      == "chars {B, A, C, D, E, F}, links {[0, 1], [2, 3], [4, 5]}"
    );

    instanceCopy.columnSwap(3, 4);
    BOOST_CHECK(
      instanceCopy.toString()
      == "chars {B, A, C, E, D, F}, links {[0, 1], [2, 4], [3, 5]}"
    );

    instanceCopy.columnSwap(0, 5);
    BOOST_CHECK(
      instanceCopy.toString()
      == "chars {F, A, C, E, D, B}, links {[0, 3], [1, 5], [2, 4]}"
    );
  }

  { // columnSmaller
    auto instanceCopy = instanceWithBondedLigands;
    BOOST_CHECK(
      instanceCopy.columnSmaller(0, 1)
      && instanceCopy.columnSmaller(1, 2)
      && instanceCopy.columnSmaller(2, 3)
      && instanceCopy.columnSmaller(3, 4)
      && instanceCopy.columnSmaller(4, 5)
    );
  }

  { // reverseColumns
    auto instanceCopy = instanceWithBondedLigands;
    instanceCopy.reverseColumns(0, 6);
    BOOST_CHECK(
      instanceCopy.toString()
      == "chars {F, E, D, C, B, A}, links {[0, 1], [2, 3], [4, 5]}"
    );

    instanceCopy.reverseColumns(2, 6);
    BOOST_CHECK(
      instanceCopy.toString()
      == "chars {F, E, A, B, C, D}, links {[0, 1], [2, 3], [4, 5]}"
    );
  }
}

BOOST_AUTO_TEST_CASE( columnSmallerConsistency ) {
  Stereopermutation single {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'A', 'A', 'A', 'A'},
    {
      std::make_pair(0,1),
      std::make_pair(2,3),
      std::make_pair(4,5)
    }
  };

  bool pass = true;

  do {
    for(unsigned i = 0; i < 6 && pass; i++) {
      for(unsigned j = i + 1; j < 6 && pass; j++) {
        if(!OperatorTests::XOR(
            single.columnSmaller(i, j),
            single.columnSmaller(j, i),
            !single.columnSmaller(i, j) && !single.columnSmaller(j, i)
        )) {
          pass = false;
        }
      }
    }
  } while(single.nextPermutation() && pass);

  BOOST_CHECK(pass);
}

BOOST_AUTO_TEST_CASE( rotationCorrectness ) {
  Stereopermutation testCase {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'C', 'D', 'B', 'B'},
    {
      std::make_pair(0, 5),
      std::make_pair(1, 4)
    }
  };

  auto isAorB = [](const char& test) -> bool {
    return (
      test == 'A'
      || test == 'B'
    );
  };

  auto testInstance = [&isAorB](
    const Stereopermutation& instance
  ) {
    return std::accumulate(
      instance.links.begin(),
      instance.links.end(),
      true,
      [&isAorB, &instance](
        const bool& carry,
        const std::pair<unsigned, unsigned>& pair
      ) {
        return (
          carry
          && isAorB(
            instance.characters.at(
              pair.first
            )
          ) && isAorB(
            instance.characters.at(
              pair.second
            )
          )
        );
      }
    );
  };

  auto allRotations = testCase.generateAllRotations(Symmetry::Name::Octahedral);
  for(const auto& copy : allRotations) {
    BOOST_CHECK(testInstance(copy));
  }
}

BOOST_AUTO_TEST_CASE( octahedralSymmetryCorrectness ) {
  Stereopermutation octahedralInstance(
    Symmetry::Name::Octahedral,
    {'A', 'B', 'C', 'D', 'E', 'F'}
  );

  BOOST_CHECK(
    octahedralInstance.generateAllRotations(Symmetry::Name::Octahedral).size()
    == 24 // 4 C4 cases on each of 6 A position selections
  );

}

void run_tests_with_counts(
  const Symmetry::Name& symmetryName,
  const std::vector<
    std::tuple<
      std::vector<char>, // characters
      std::set<
        std::pair<unsigned, unsigned>
      >, // pairs
      unsigned // expectedUnique
    >
  >& test_cases
) {
  for(const auto& tuple: test_cases) {
    std::vector<char> characters;
    std::set<
      std::pair<unsigned, unsigned>
    > pairs;
    unsigned expectedUnique;

    std::tie(characters, pairs, expectedUnique) = tuple;

    // instantiate
    Stereopermutation assignment = (pairs.size() == 0)
      ? Stereopermutation(symmetryName, characters)
      : Stereopermutation(symmetryName, characters, pairs);

    auto unique = uniquesWithWeights(assignment, symmetryName);

    BOOST_CHECK(unique.assignments.size() == expectedUnique );

    if(unique.assignments.size() != expectedUnique) {
      std::cout << "Mismatch: Expected " << expectedUnique
        << " assignments for: \n" << assignment << ", got "
        << unique.assignments.size() << " assignments:" << std::endl;
    }

    std::cout << "{";
    for(unsigned i = 0; i < characters.size(); i++) {
      std::cout << characters[i];
      if(i != characters.size() - 1) std::cout << ", ";
    }
    std::cout << "} " << Symmetry::name(symmetryName) << std::endl;

    for(unsigned i = 0; i < unique.assignments.size(); i++) {
      std::cout << "Weight " << unique.weights[i] << ": "
        << unique.assignments[i] << std::endl;
    }

    std::cout << "----------------------------------------" << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( individual_bugfixes ) {
  Stereopermutation a {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(2, 3),
      std::make_pair(1, 4),
      std::make_pair(0, 5)
    }
  };
  Stereopermutation b {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'B', 'A', 'B', 'B'},
    {
      std::make_pair(3, 5),
      std::make_pair(1, 2),
      std::make_pair(0, 4)
    }
  };

  // one and only one of the following can be true for any Stereopermutations a and b
  BOOST_CHECK(OperatorTests::testLogicalOperators(a, b));

  BOOST_CHECK(a.isRotationallySuperimposable(b, Symmetry::Name::Octahedral));
  BOOST_CHECK(b.isRotationallySuperimposable(a, Symmetry::Name::Octahedral));

  /* Contrived example of two that have inconsistent logical operators, just
   * reordered op pairs. Will evaluate == but also < w/ current impl.
   */
  Stereopermutation c {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(2, 3),
      std::make_pair(1, 4),
      std::make_pair(0, 5)
    }
  };
  Stereopermutation d {
    Symmetry::Name::Octahedral,
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(1, 4),
      std::make_pair(2, 3),
      std::make_pair(0, 5)
    }
  };

  BOOST_CHECK(OperatorTests::testLogicalOperators(c, d));
  BOOST_CHECK(c.isRotationallySuperimposable(d, Symmetry::Name::Octahedral));
  BOOST_CHECK(d.isRotationallySuperimposable(c, Symmetry::Name::Octahedral));
}

/* Tetrahedral tests */
/* These were thought up myself */
BOOST_AUTO_TEST_CASE( tetrahedral_monodentate ) {
  run_tests_with_counts(
    Symmetry::Name::Tetrahedral,
    {
    // M_A
    std::make_tuple(
      std::vector<char>(4, 'A'),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      1
    ),
    // M_A3B
    std::make_tuple(
      std::vector<char>({'A', 'A', 'A', 'B'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      1
    ),
    // M_A2B2
    std::make_tuple(
      std::vector<char>({'A', 'A', 'B', 'B'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      1
    ),
    // M_A2BC
    std::make_tuple(
      std::vector<char>({'A', 'A', 'B', 'C'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      1
    ),
    // M_ABCD
    std::make_tuple(
      std::vector<char>({'A', 'B', 'C', 'D'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      2
    )
  });
}

/* Square Planar tests */
/* These were thought up myself */
BOOST_AUTO_TEST_CASE( square_planar_monodentate ) {
  run_tests_with_counts(
    Symmetry::Name::SquarePlanar,
    {
    // M_A
    std::make_tuple(
      std::vector<char>(4, 'A'),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      1
    ),
    // M_A3B
    std::make_tuple(
      std::vector<char>({'A', 'A', 'A', 'B'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      1
    ),
    // M_A2B2
    std::make_tuple(
      std::vector<char>({'A', 'A', 'B', 'B'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      2
    ),
    // M_A2BC
    std::make_tuple(
      std::vector<char>({'A', 'A', 'B', 'C'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      2
    ),
    // M_ABCD
    std::make_tuple(
      std::vector<char>({'A', 'B', 'C', 'D'}),
      std::set<
        std::pair<unsigned, unsigned>
      >(),
      3
    )
  });
}

/* Octahedral tests */
/* Expected values taken from
 * Miessler, Gary L., Tarr, Donald A.: Inorganic Chemistry, Third Edition.
 * Do not know whether these are correct! They are "all calculated using a
 * computer program [...]."
 * The reference however is useful: WE Bennett, Inorg. Chem. 1969
 */
BOOST_AUTO_TEST_CASE( octahedral_monodentate ) {
  run_tests_with_counts(
    Symmetry::Name::Octahedral,
    {
      // M_A
      std::make_tuple(
        std::vector<char>(6, 'A'),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        1
      ),
      // M_AB
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'A', 'B'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        1
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'B', 'B'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        2
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'B', 'B', 'B'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        2
      ),
      // M_ABC
      /* A4 B C
       * A3 B2 C
       * not A3 B C2 == A3 B2 C
       * A2 B2 C2
       * not A2 B3 C == A3 B2 C
       */
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'B', 'C'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        2
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'B', 'B', 'C'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        3
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'B', 'C', 'C'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        6
      ),
      // M_ABCD
      /* A3 B C D
       * A2 B2 C D
       */
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'B', 'C', 'D'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        5
      ),
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'B', 'C', 'D'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        8
      ),
      // M_ABCDE
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'C', 'D', 'E'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        15
      ),
      // M_ABCDEF
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'D', 'E', 'F'}),
        std::set<
          std::pair<unsigned, unsigned>
        >(),
        30
      )
    }
  );
}

BOOST_AUTO_TEST_CASE( octahedral_multidentate ) {
  run_tests_with_counts(
    Symmetry::Name::Octahedral,
    {
      // M(A-A)_3
      std::make_tuple(
        std::vector<char>({'A', 'A', 'A', 'A', 'A', 'A'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3),
          std::make_pair(4, 5)
        }),
        2
        /* NOTE: besides the two cis-cis-cis enantiomers, there are
         * cis-cis-trans and trans-trans-trans isomers, these are removed by
         * default!
         */
      ),
      // M(A-B)_3
      std::make_tuple(
        std::vector<char>({'A', 'B', 'A', 'B', 'A', 'B'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3),
          std::make_pair(4, 5)
        }),
        4
      ),
      // M(A-B)_2 CD
      std::make_tuple(
        std::vector<char>({'A', 'B', 'A', 'B', 'C', 'D'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3)
        }),
        11
      ),
      // M(A-A)(B-C)DE
      std::make_tuple(
        std::vector<char>({'A', 'A', 'B', 'C', 'D', 'E'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3)
        }),
        10
      ),
      // M(A-B)(C-D)EF
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'D', 'E', 'F'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(2, 3)
        }),
        20
      ),
      // M(A-B-A)CDE
      std::make_tuple(
        std::vector<char>({'A', 'B', 'A', 'C', 'D', 'E'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2)
        }),
        9
      ),
      // M(A-B-C)_2
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'A', 'B', 'C'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2),
          std::make_pair(3, 4),
          std::make_pair(4, 5)
        }),
        11
      ),
      // M(A-B-B-A)CD
      std::make_tuple(
        std::vector<char>({'A', 'B', 'B', 'A', 'C', 'D'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2),
          std::make_pair(2, 3)
        }),
        7
      ),
      // M(A-B-C-B-A)D
      std::make_tuple(
        std::vector<char>({'A', 'B', 'C', 'B', 'A', 'D'}),
        std::set<
          std::pair<unsigned, unsigned>
        >({
          std::make_pair(0, 1),
          std::make_pair(1, 2),
          std::make_pair(2, 3),
          std::make_pair(3, 4)
        }),
        7
      )
    }
  );
}
