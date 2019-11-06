/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#define BOOST_TEST_MODULE ConnectivityManagerTests
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <iostream>
#include <cassert>
#include <functional>
#include <numeric>

#include "shapes/Data.h"
#include "shapes/Properties.h"

#include "stereopermutation/GenerateUniques.h"
#include "stereopermutation/Composites.h"

#include "temple/Functional.h"
#include "temple/Random.h"
#include "temple/Stringify.h"
#include "temple/constexpr/LogicalOperatorTests.h"

using namespace Scine;
using namespace stereopermutation;

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
BOOST_AUTO_TEST_CASE(StereopermutationInstantiation) {
  for(const auto& shape: Shapes::allShapes) {
    const unsigned S = Shapes::size(shape);
    if(S > 8) {
      continue;
    }

    Stereopermutation testStereopermutation(std::vector<char>(S, 'A'));
  }
}

BOOST_AUTO_TEST_CASE(StereopermutationBasics) {
  // Constructors
  Stereopermutation instanceWithBondedLigands(
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
    const auto& instanceCopy = instanceWithBondedLigands;
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

BOOST_AUTO_TEST_CASE(ColumnComparisonConsistency) {
  Stereopermutation single {
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
        if(!temple::Math::XOR(
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
        const bool carry,
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

  auto allRotations = testCase.generateAllRotations(Shapes::Shape::Octahedron);
  for(const auto& copy : allRotations) {
    BOOST_CHECK(testInstance(copy));
  }
}

BOOST_AUTO_TEST_CASE(OctahedralSymmetryCorrectness) {
  Stereopermutation octahedralInstance(
    {'A', 'B', 'C', 'D', 'E', 'F'}
  );

  BOOST_CHECK(
    octahedralInstance.generateAllRotations(Shapes::Shape::Octahedron).size()
    == 24 // 4 C4 cases on each of 6 A position selections
  );

}

void run_tests_with_counts(
  const Shapes::Shape& shape,
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
    Stereopermutation stereopermutation = pairs.empty()
      ? Stereopermutation(characters)
      : Stereopermutation(characters, pairs);

    // The count of uniques in the tests are all without trans-arranged pairs!
    auto unique = uniquesWithWeights(
      stereopermutation,
      shape,
      true // remove trans arranged pairs
    );

    BOOST_CHECK(unique.stereopermutations.size() == expectedUnique );

    if(unique.stereopermutations.size() != expectedUnique) {
      std::cout << "Mismatch: Expected " << expectedUnique
        << " stereopermutations for: \n" << stereopermutation << ", got "
        << unique.stereopermutations.size() << " stereopermutations:" << std::endl;
    }
  }
}

BOOST_AUTO_TEST_CASE(BugfixTests) {
  Stereopermutation a {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(2, 3),
      std::make_pair(1, 4),
      std::make_pair(0, 5)
    }
  };
  Stereopermutation b {
    {'A', 'A', 'B', 'A', 'B', 'B'},
    {
      std::make_pair(3, 5),
      std::make_pair(1, 2),
      std::make_pair(0, 4)
    }
  };

  // one and only one of the following can be true for any Stereopermutations a and b
  BOOST_CHECK(temple::testLogicalOperators(a, b));

  BOOST_CHECK(a.isRotationallySuperimposable(b, Shapes::Shape::Octahedron));
  BOOST_CHECK(b.isRotationallySuperimposable(a, Shapes::Shape::Octahedron));

  /* Contrived example of two that have inconsistent logical operators, just
   * reordered op pairs. Will evaluate == but also < w/ current impl.
   */
  Stereopermutation c {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(2, 3),
      std::make_pair(1, 4),
      std::make_pair(0, 5)
    }
  };
  Stereopermutation d {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {
      std::make_pair(1, 4),
      std::make_pair(2, 3),
      std::make_pair(0, 5)
    }
  };

  BOOST_CHECK(temple::testLogicalOperators(c, d));
  BOOST_CHECK(c.isRotationallySuperimposable(d, Shapes::Shape::Octahedron));
  BOOST_CHECK(d.isRotationallySuperimposable(c, Shapes::Shape::Octahedron));
}

/* Tetrahedron tests */
BOOST_AUTO_TEST_CASE(MonodentateTetrahedral) {
  run_tests_with_counts(
    Shapes::Shape::Tetrahedron,
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
BOOST_AUTO_TEST_CASE(MonodentateSquarePlanar) {
  run_tests_with_counts(
    Shapes::Shape::Square,
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

/* Octahedron tests */
/* Expected values taken from
 * Miessler, Gary L., Tarr, Donald A.: Inorganic Chemistry, Third Edition.
 * Do not know whether these are correct! They are "all calculated using a
 * computer program [...]."
 * The reference however is useful: WE Bennett, Inorg. Chem. 1969
 */
BOOST_AUTO_TEST_CASE(MonodentateOctahedron) {
  run_tests_with_counts(
    Shapes::Shape::Octahedron,
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

BOOST_AUTO_TEST_CASE(MultidentateOctahedron) {
  run_tests_with_counts(
    Shapes::Shape::Octahedron,
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

void writeState(const Composite::OrientationState& a) {
  std::cout << "OrientationState {\n"
    << "  symmetry: " << Shapes::name(a.shape) << "\n"
    << "  fusedPosition: " << a.fusedVertex << "\n"
    << "  characters: " << temple::stringify(a.characters) << "\n}\n";
}

bool testOrientationState(const Composite::OrientationState& a) {
  // Make a copy and modify that
  auto aCopy = a;

  // Transform and revert the OrientationState
  auto reversionMapping = aCopy.transformToCanonical();
  aCopy.revert(reversionMapping);

  return aCopy == a;
}

// Ensure that transformation and reversion work the way they should
BOOST_AUTO_TEST_CASE(OrientationStateTests) {
  for(const auto& shape : Shapes::allShapes) {
    const unsigned S = Shapes::size(shape);
    if(S > 8) {
      continue;
    }

    std::vector<char> maximumAsymmetricCase (S);
    for(unsigned i = 0; i < S; ++i) {
      maximumAsymmetricCase.at(i) = 'A' + i;
    }

    for(unsigned i = 0; i < S; ++i ) {
      BOOST_CHECK_MESSAGE(
        testOrientationState(
          Composite::OrientationState {
            shape,
            i,
            maximumAsymmetricCase,
            0
          }
        ),
        "Transformation and reversion does not work for "
        << Shapes::name(shape)
        << " on position " << i << "."
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(compositesTests) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  Composite a {
    Composite::OrientationState {
      Shapes::Shape::Seesaw,
      0,
      {'A', 'B', 'C', 'D'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0,
      {'A', 'B', 'C', 'D'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    a.permutations() == 3u,
    "Expected 3 permutations, got " << a.permutations()
  );

  Composite b {
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      4,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      2,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    b.permutations() == 4u,
    "Expected 4 permutations, got " << b.permutations()
  );

  Composite c {
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::EquilateralTriangle,
      0,
      {'A', 'B', 'C'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    c.permutations() == 2,
    "Expected 2 permutations, got " << c.permutations()
  );

  Composite d {
    Composite::OrientationState {
      Shapes::Shape::EquilateralTriangle,
      0,
      {'A', 'B', 'C'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::EquilateralTriangle,
      0,
      {'A', 'B', 'C'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    d.permutations() == 2,
    "Expected 2 permutations, got " << d.permutations()
  );

  Composite e {
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0,
      {'A', 'B'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    e.permutations() == 2,
    "Expected 2 permutations, got " << e.permutations()
  );
}

BOOST_AUTO_TEST_CASE(compositesAlignment) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  Composite a {
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0,
      {'A', 'B', 'C', 'D'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0,
      {'A', 'B', 'C', 'D'},
      rightIdentifier
    },
    Composite::Alignment::Staggered
  };

  BOOST_CHECK_MESSAGE(
    a.permutations() == 3,
    "Expected three permutations, got " << a.permutations()
  );

  Composite b {
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::SquarePyramid,
      0,
      {'A', 'B', 'C', 'D', 'E'},
      rightIdentifier
    },
    Composite::Alignment::Staggered
  };

  BOOST_CHECK_MESSAGE(
    b.permutations() == 3,
    "Expected three permutations, got " << b.permutations()
  );

  Composite c {
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      4,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      2,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    c.permutations() == 4,
    "Expected four permutations, got " << c.permutations()
  );
}

BOOST_AUTO_TEST_CASE(numUnlinkedStereopermutationsTest) {
  // Crosscheck number of unlinked stereopermutations with shapes
  for(const Shapes::Shape shape : Shapes::allShapes) {
    const unsigned S = Shapes::size(shape);
    if(S > 8) {
      continue;
    }

    std::vector<char> characters (S);
    for(unsigned nIdentical = 1; nIdentical < S; ++nIdentical) {
      // Populate characters for construction of a Stereopermutation
      for(unsigned i = 0; i < characters.size(); ++i) {
        characters.at(i) = 'A' + std::max(
          0,
          static_cast<int>(i) - static_cast<int>(nIdentical) + 1
        );
      }

      assert(
        std::count(
          std::begin(characters),
          std::end(characters),
          'A'
        ) == nIdentical
      );

      const Stereopermutation initialStereopermutation {characters};

      const unsigned uniquesCount = uniques(initialStereopermutation, shape).size();

      // shapes prediction
      const unsigned chemicalSymmetriesCount = Shapes::properties::numUnlinkedStereopermutations(
        shape,
        nIdentical
      );

      BOOST_CHECK_MESSAGE(
        uniquesCount == chemicalSymmetriesCount,
        "stereopermutation and shapes differ in unique "
        << "stereopermutation counts for " << Shapes::name(shape)
        << " and " << nIdentical << " ligands: " << uniquesCount << " vs "
        << chemicalSymmetriesCount
      );
    }
  }
}
