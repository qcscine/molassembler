/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <cassert>
#include <functional>
#include <numeric>

#include "molassembler/Shapes/Data.h"
#include "molassembler/Shapes/Properties.h"

#include "molassembler/Stereopermutation/Manipulation.h"
#include "molassembler/Stereopermutation/Composites.h"
#include "molassembler/Stereopermutation/RotationEnumerator.h"

#include "molassembler/Temple/Functional.h"
#include "molassembler/Temple/Random.h"
#include "molassembler/Temple/Stringify.h"
#include "molassembler/Temple/constexpr/LogicalOperatorTests.h"

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

shapes::Vertex operator "" _v(unsigned long long v) {
  return shapes::Vertex(v);
}

// create instances of all symmetries with monodentate ligands
BOOST_AUTO_TEST_CASE(StereopermutationInstantiation) {
  for(const auto& shape: shapes::allShapes) {
    const unsigned S = shapes::size(shape);
    if(S > 8) {
      continue;
    }

    Stereopermutation testStereopermutation(std::vector<char>(S, 'A'));
  }
}

BOOST_AUTO_TEST_CASE(StereopermutationRotation) {
  // Octahedron ship-screw like cis-cis-cis
  const Stereopermutation asymm {
    {'A', 'B', 'C', 'D', 'E', 'F'},
    {{0_v, 1_v}, {2_v, 4_v}, {3_v, 5_v}}
  };

  std::vector<shapes::Vertex> C4z {{3_v, 0_v, 1_v, 2_v, 4_v, 5_v}};

  const Stereopermutation expected {
    {'D', 'A', 'B', 'C', 'E', 'F'},
    {{0_v, 5_v}, {1_v, 2_v}, {3_v, 4_v}}
  };

  const Stereopermutation rotated = asymm.applyPermutation(C4z);

  BOOST_CHECK_MESSAGE(
    rotated == expected,
    "Expected rotation of " << asymm.toString() << " with C4z in octahedron to\n"
    << expected.toString() << ", got\n" << rotated.toString()
  );
}

BOOST_AUTO_TEST_CASE(RotationCorrectness) {
  Stereopermutation testCase {
    {'A', 'A', 'C', 'D', 'B', 'B'},
    {{0_v, 5_v}, {1_v, 4_v}}
  };

  auto isAorB = [](const char& test) -> bool {
    return (test == 'A' || test == 'B');
  };

  const auto allRotations = generateAllRotations(testCase, shapes::Shape::Octahedron);
  // Make sure the links still refer to A or B characters
  for(const auto& rotation : allRotations) {
    BOOST_CHECK(
      temple::all_of(
        rotation.links,
        [&](const auto& link) -> bool {
          return (
            isAorB(rotation.characters.at(link.first))
            && isAorB(rotation.characters.at(link.second))
          );
        }
      )
    );
  }
}

BOOST_AUTO_TEST_CASE(OctahedralSymmetryCorrectness) {
  Stereopermutation octahedralInstance(
    {'A', 'B', 'C', 'D', 'E', 'F'}
  );

  BOOST_CHECK_EQUAL(
    generateAllRotations(octahedralInstance, shapes::Shape::Octahedron).size(),
    24 // 4 C4 cases on each of 6 A position selections
  );
}

BOOST_AUTO_TEST_CASE(BugfixTests) {
  Stereopermutation a {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {{2_v, 3_v}, {1_v, 4_v}, {0_v, 5_v}}
  };
  Stereopermutation b {
    {'A', 'A', 'B', 'A', 'B', 'B'},
    {{3_v, 5_v}, {1_v, 2_v}, {0_v, 4_v}}
  };

  // one and only one of the following can be true for any Stereopermutations a and b
  BOOST_CHECK(temple::testLogicalOperators(a, b));

  BOOST_CHECK(rotationallySuperimposable(a, b, shapes::Shape::Octahedron));
  BOOST_CHECK(rotationallySuperimposable(b, a, shapes::Shape::Octahedron));

  /* Contrived example of two that have inconsistent logical operators, just
   * reordered op pairs. Will evaluate == but also < w/ current impl.
   */
  Stereopermutation c {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {{2_v, 3_v}, {1_v, 4_v}, {0_v, 5_v}}
  };
  Stereopermutation d {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {{1_v, 4_v}, {2_v, 3_v}, {0_v, 5_v}}
  };

  BOOST_CHECK(temple::testLogicalOperators(c, d));
  BOOST_CHECK(rotationallySuperimposable(c, d, shapes::Shape::Octahedron));
  BOOST_CHECK(rotationallySuperimposable(d, c, shapes::Shape::Octahedron));
}

using Characters = Stereopermutation::CharacterOccupation;
using PairSet = std::vector<
  std::pair<unsigned, unsigned>
>;

struct TestCase {
  Characters characters;
  PairSet pairs;
  unsigned expectedUnique;
};

Stereopermutation::OrderedLinks makeLinks(const PairSet& pairs) {
  return temple::map(
    pairs,
    [](const std::pair<unsigned, unsigned>& v) -> Stereopermutation::Link {
      return std::make_pair(
        shapes::Vertex(v.first),
        shapes::Vertex(v.second)
      );
    }
  );
}

void runTestsWithCounts(
  const shapes::Shape& shape,
  const std::vector<TestCase>& testCases
) {
  for(const auto& testCase: testCases) {
    // instantiate
    const Stereopermutation stereopermutation {
      testCase.characters,
      makeLinks(testCase.pairs)
    };

    // The count of uniques in the tests are all without trans-arranged pairs!
    const auto unique = uniques(
      stereopermutation,
      shape,
      true // remove trans arranged pairs
    );

    // std::cout << stereopermutation.toString() << " unique weights: " << temple::stringify(unique.weights) << "\n";

    BOOST_CHECK_MESSAGE(
      unique.list.size() == testCase.expectedUnique,
      "Mismatch: Expected " << testCase.expectedUnique
        << " stereopermutations for\n" << stereopermutation.toString() << " in " << shapes::name(shape) << " shape, got "
        << unique.list.size() << " stereopermutations\n"
    );
  }
}

/* Tetrahedron tests */
BOOST_AUTO_TEST_CASE(MonodentateTetrahedral) {
  runTestsWithCounts(
    shapes::Shape::Tetrahedron,
    {
      { // AAAA
        Characters(4, 'A'),
        {},
        1
      },
      { // A3B
        Characters({'A', 'A', 'A', 'B'}),
        {},
        1
      },
      { // A2B2
        Characters({'A', 'A', 'B', 'B'}),
        {},
        1
      },
      { // A2BC
        Characters({'A', 'A', 'B', 'C'}),
        {},
        1
      },
      { // ABCD
        Characters({'A', 'B', 'C', 'D'}),
        {},
        2
      },
    }
  );
}

/* Square Planar tests */
BOOST_AUTO_TEST_CASE(MonodentateSquarePlanar) {
  runTestsWithCounts(
    shapes::Shape::Square,
    {
      { // M_A
        Characters(4, 'A'),
        {},
        1
      },
      { // M_A3B
        Characters({'A', 'A', 'A', 'B'}),
        {},
        1
      },
      { // M_A2B2
        Characters({'A', 'A', 'B', 'B'}),
        {},
        2
      },
      { // M_A2BC
        Characters({'A', 'A', 'B', 'C'}),
        {},
        2
      },
      { // M_ABCD
        Characters({'A', 'B', 'C', 'D'}),
        {},
        3
      },
    }
  );
}

/* Octahedron tests */
/* Expected values taken from
 * Miessler, Gary L., Tarr, Donald A.: Inorganic Chemistry, Third Edition.
 * Do not know whether these are correct! They are "all calculated using a
 * computer program [...]."
 * The reference however is useful: WE Bennett, Inorg. Chem. 1969
 */
BOOST_AUTO_TEST_CASE(MonodentateOctahedron) {
  runTestsWithCounts(
    shapes::Shape::Octahedron,
    {
      {
        Characters(6, 'A'),
        {},
        1
      },
      {
        Characters({'A', 'A', 'A', 'A', 'A', 'B'}),
        {},
        1
      },
      {
        Characters({'A', 'A', 'A', 'A', 'B', 'B'}),
        {},
        2
      },
      {
        Characters({'A', 'A', 'A', 'B', 'B', 'B'}),
        {},
        2
      },
      {
        Characters({'A', 'A', 'A', 'A', 'B', 'C'}),
        {},
        2
      },
      {
        Characters({'A', 'A', 'A', 'B', 'B', 'C'}),
        {},
        3
      },
      {
        Characters({'A', 'A', 'B', 'B', 'C', 'C'}),
        {},
        6
      },
      {
        Characters({'A', 'A', 'A', 'B', 'C', 'D'}),
        {},
        5
      },
      {
        Characters({'A', 'A', 'B', 'B', 'C', 'D'}),
        {},
        8
      },
      {
        Characters({'A', 'A', 'B', 'C', 'D', 'E'}),
        {},
        15
      },
      {
        Characters({'A', 'B', 'C', 'D', 'E', 'F'}),
        {},
        30
      }
    }
  );
}

BOOST_AUTO_TEST_CASE(MultidentateOctahedron) {
  runTestsWithCounts(
    shapes::Shape::Octahedron,
    {
      { // M(A-A)_3
        Characters({'A', 'A', 'A', 'A', 'A', 'A'}),
        PairSet({{0, 1}, {2, 3}, {4, 5}}),
        2
        /* NOTE: besides the two cis-cis-cis enantiomers, there are
         * cis-cis-trans and trans-trans-trans isomers, these are removed by
         * default!
         */
      },
      { // M(A-B)_3
        Characters({'A', 'B', 'A', 'B', 'A', 'B'}),
        PairSet({{0, 1}, {2, 3}, {4, 5}}),
        4
      },
      { // M(A-B)_2 CD
        Characters({'A', 'B', 'A', 'B', 'C', 'D'}),
        PairSet({{0, 1}, {2, 3}}),
        11
      },
      { // M(A-A)(B-C)DE
        Characters({'A', 'A', 'B', 'C', 'D', 'E'}),
        PairSet({{0, 1}, {2, 3}}),
        10
      },
      { // M(A-B)(C-D)EF
        Characters({'A', 'B', 'C', 'D', 'E', 'F'}),
        PairSet({{0, 1}, {2, 3}}),
        20
      },
      { // M(A-B-A)CDE
        Characters({'A', 'B', 'A', 'C', 'D', 'E'}),
        PairSet({{0, 1}, {1, 2}}),
        9
      },
      { // M(A-B-C)_2
        Characters({'A', 'B', 'C', 'A', 'B', 'C'}),
        PairSet({{0, 1}, {1, 2}, {3, 4}, {4, 5}}),
        11
      },
      { // M(A-B-B-A)CD
        Characters({'A', 'B', 'B', 'A', 'C', 'D'}),
        PairSet({{0, 1}, {1, 2}, {2, 3}}),
        7
      },
      { // M(A-B-C-B-A)D
        Characters({'A', 'B', 'C', 'B', 'A', 'D'}),
        PairSet({{0, 1}, {1, 2}, {2, 3}, {3, 4}}),
        7
      }
    }
  );
}

BOOST_AUTO_TEST_CASE(OctahedralCaseWithoutTransRemoval) {
  const auto shape = shapes::Shape::Octahedron;
  const auto characters = Characters(6, 'A');
  const auto pairs = PairSet({{0, 1}, {2, 3}, {4, 5}});

  const Stereopermutation stereopermutation {characters, makeLinks(pairs)};
  const auto unique = uniques(stereopermutation, shape);
  BOOST_CHECK_EQUAL(unique.list.size(), 4u);
}

bool testOrientationState(Composite::OrientationState a) {
  // Make a copy and modify that
  auto aCopy = a;

  // Transform and revert the OrientationState
  auto reversionMapping = aCopy.transformToCanonical();
  aCopy.revert(reversionMapping);

  return aCopy == a;
}

// Ensure that transformation and reversion work the way they should
BOOST_AUTO_TEST_CASE(OrientationStateTests) {
  for(const auto& shape : shapes::allShapes) {
    const unsigned S = shapes::size(shape);
    if(S > 8) {
      continue;
    }

    std::vector<char> maximumAsymmetricCase (S);
    for(unsigned i = 0; i < S; ++i) {
      maximumAsymmetricCase.at(i) = 'A' + i;
    }

    for(shapes::Vertex i {0}; i < S; ++i ) {
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
        << shapes::name(shape)
        << " on position " << i << "."
      );
    }
  }
}

BOOST_AUTO_TEST_CASE(CompositeExamples) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  Composite a {
    Composite::OrientationState {
      shapes::Shape::Seesaw,
      0_v,
      {'A', 'B', 'C', 'D'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::Tetrahedron,
      0_v,
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
      shapes::Shape::Octahedron,
      4_v,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::Octahedron,
      2_v,
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
      shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::EquilateralTriangle,
      0_v,
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
      shapes::Shape::EquilateralTriangle,
      0_v,
      {'A', 'B', 'C'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::EquilateralTriangle,
      0_v,
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
      shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    e.permutations() == 2,
    "Expected 2 permutations, got " << e.permutations()
  );

  Composite f {
    Composite::OrientationState {
      shapes::Shape::Line,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::Line,
      0_v,
      {'A', 'B'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    f.permutations() == 0,
    "Expected 0 permutations, got " << f.permutations()
  );
}

BOOST_AUTO_TEST_CASE(CompositeAlignment) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  Composite a {
    Composite::OrientationState {
      shapes::Shape::Tetrahedron,
      0_v,
      {'A', 'B', 'C', 'D'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::Tetrahedron,
      0_v,
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
      shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::SquarePyramid,
      0_v,
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
      shapes::Shape::Octahedron,
      4_v,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      leftIdentifier
    },
    Composite::OrientationState {
      shapes::Shape::Octahedron,
      2_v,
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
  for(const shapes::Shape shape : shapes::allShapes) {
    const unsigned S = shapes::size(shape);
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

      const unsigned uniquesCount = uniques(initialStereopermutation, shape).list.size();

      // shapes prediction
      const unsigned chemicalSymmetriesCount = shapes::properties::numUnlinkedStereopermutations(
        shape,
        nIdentical
      );

      BOOST_CHECK_MESSAGE(
        uniquesCount == chemicalSymmetriesCount,
        "stereopermutation and shapes differ in unique "
        << "stereopermutation counts for " << shapes::name(shape)
        << " and " << nIdentical << " ligands: " << uniquesCount << " vs "
        << chemicalSymmetriesCount
      );
    }
  }
}
