#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>

#include "GenerateUniques.h"
#include "SymmetryInformation.h"

#include "LogicalOperatorTests.h"

using namespace UniqueAssignments;

/* TODO
 * - investigate failures!
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

BOOST_AUTO_TEST_CASE( assignment_instantiation ) {
  // can make instances of all symmetries
  Assignment<PermSymmetry::Tetrahedral> tetr(
    std::vector<char>(4, 'A')
  );
  Assignment<PermSymmetry::SquarePlanar> sqpl(
    std::vector<char>(4, 'A')
  );
  Assignment<PermSymmetry::SquarePyramidal> sqpy(
    std::vector<char>(5, 'A')
  );
  Assignment<PermSymmetry::TrigonalBiPyramidal> trigbipy(
    std::vector<char>(5, 'A')
  );
  Assignment<PermSymmetry::Octahedral> octa(
    std::vector<char>(6, 'A')
  );
}

BOOST_AUTO_TEST_CASE( assignment_basics ) {
  // Constructors
  Assignment<PermSymmetry::Octahedral> instanceWithBondedLigands(
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

BOOST_AUTO_TEST_CASE( octahedralSymmetryCorrectness ) {
  Assignment<PermSymmetry::Octahedral> octahedralInstance(
    {'A', 'B', 'C', 'D', 'E', 'F'}
  );

  BOOST_CHECK(
    octahedralInstance.generateAllRotations().size()
    == 24 // 4 C4 cases on each of 6 A position selections
  );
  
}

template<class Symmetry>
void run_tests(
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
    Assignment<Symmetry> assignment = (pairs.size() == 0)
      ? Assignment<Symmetry>(characters)
      : Assignment<Symmetry>(characters, pairs);

    auto unique = uniqueAssignments(assignment);
    BOOST_CHECK(unique.size() == expectedUnique );
    if(unique.size() != expectedUnique) {
      std::cout << "Mismatch: Expected " << expectedUnique
        << " assignments for: \n" << assignment << ", got " 
        << unique.size() << " assignments:" << std::endl;
      for(const auto& uniqueAssignment: unique) {
        std::cout << uniqueAssignment << std::endl;
      }
    } 
  }
}

BOOST_AUTO_TEST_CASE( individual_bugfixes ) {
  Assignment<PermSymmetry::Octahedral> a {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(2, 3),
      std::make_pair(1, 4),
      std::make_pair(0, 5)
    }
  };
  Assignment<PermSymmetry::Octahedral> b {
    {'A', 'A', 'B', 'A', 'B', 'B'},
    {
      std::make_pair(3, 5),
      std::make_pair(1, 2),
      std::make_pair(0, 4)
    }
  };

  // one and only one of the following can be true for any Assignments a and b
  BOOST_CHECK(OperatorTests::testLogicalOperators(a, b));

  BOOST_CHECK(a.isRotationallySuperimposable(b));
  BOOST_CHECK(b.isRotationallySuperimposable(a));

  /* Contrived example of two that have inconsistent logical operators, just
   * reordered op pairs. Will evaluate == but also < w/ current impl.
   */
  Assignment<PermSymmetry::Octahedral> c {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(2, 3),
      std::make_pair(1, 4),
      std::make_pair(0, 5)
    }
  };
  Assignment<PermSymmetry::Octahedral> d {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(1, 4),
      std::make_pair(2, 3),
      std::make_pair(0, 5)
    }
  };

  BOOST_CHECK(OperatorTests::testLogicalOperators(c, d));
  BOOST_CHECK(c.isRotationallySuperimposable(d));
  BOOST_CHECK(d.isRotationallySuperimposable(c));
}

/* Tetrahedral tests */
/* These were thought up myself */
BOOST_AUTO_TEST_CASE( tetrahedral_monodentate ) {
  run_tests<PermSymmetry::Tetrahedral>({
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
  run_tests<PermSymmetry::SquarePlanar>({
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
  run_tests<PermSymmetry::Octahedral>(
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
  run_tests<PermSymmetry::Octahedral>(
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
        11 // TODO ERROR: get 10
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
        11 // TODO ERROR: get 9
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
