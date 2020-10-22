/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <vector>
#include <cassert>
#include <functional>
#include <numeric>

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Shapes/Properties.h"

#include "Molassembler/Stereopermutation/Manipulation.h"
#include "Molassembler/Stereopermutation/RotationEnumerator.h"

#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Random.h"
#include "Molassembler/Temple/Stringify.h"
#include "Molassembler/Temple/constexpr/LogicalOperatorTests.h"

using namespace Scine::Molassembler;
using namespace Stereopermutations;

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

Shapes::Vertex operator "" _v(unsigned long long v) {
  return Shapes::Vertex(v);
}

// create instances of all symmetries with monodentate ligands
BOOST_AUTO_TEST_CASE(StereopermutationInstantiation, *boost::unit_test::label("Stereopermutations")) {
  for(const auto& shape: Shapes::allShapes) {
    const unsigned S = Shapes::size(shape);
    if(S > 8) {
      continue;
    }

    Stereopermutation testStereopermutation(std::vector<char>(S, 'A'));
  }
}

BOOST_AUTO_TEST_CASE(StereopermutationRotation, *boost::unit_test::label("Stereopermutations")) {
  // Octahedron ship-screw like cis-cis-cis
  const Stereopermutation asymm {
    {'A', 'B', 'C', 'D', 'E', 'F'},
    {{0_v, 1_v}, {2_v, 4_v}, {3_v, 5_v}}
  };

  std::vector<Shapes::Vertex> C4z {{3_v, 0_v, 1_v, 2_v, 4_v, 5_v}};

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

BOOST_AUTO_TEST_CASE(RotationCorrectness, *boost::unit_test::label("Stereopermutations")) {
  Stereopermutation testCase {
    {'A', 'A', 'C', 'D', 'B', 'B'},
    {{0_v, 5_v}, {1_v, 4_v}}
  };

  auto isAorB = [](const char& test) -> bool {
    return (test == 'A' || test == 'B');
  };

  const auto allRotations = generateAllRotations(testCase, Shapes::Shape::Octahedron);
  // Make sure the links still refer to A or B characters
  for(const auto& rotation : allRotations) {
    BOOST_CHECK(
      Temple::all_of(
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

BOOST_AUTO_TEST_CASE(OctahedralSymmetryCorrectness, *boost::unit_test::label("Stereopermutations")) {
  Stereopermutation octahedralInstance(
    {'A', 'B', 'C', 'D', 'E', 'F'}
  );

  BOOST_CHECK_EQUAL(
    generateAllRotations(octahedralInstance, Shapes::Shape::Octahedron).size(),
    24 // 4 C4 cases on each of 6 A position selections
  );
}

BOOST_AUTO_TEST_CASE(BugfixTests, *boost::unit_test::label("Stereopermutations")) {
  Stereopermutation a {
    {'A', 'A', 'A', 'B', 'B', 'B'},
    {{2_v, 3_v}, {1_v, 4_v}, {0_v, 5_v}}
  };
  Stereopermutation b {
    {'A', 'A', 'B', 'A', 'B', 'B'},
    {{3_v, 5_v}, {1_v, 2_v}, {0_v, 4_v}}
  };

  // one and only one of the following can be true for any Stereopermutations a and b
  BOOST_CHECK(Temple::testLogicalOperators(a, b));

  BOOST_CHECK(rotationallySuperimposable(a, b, Shapes::Shape::Octahedron));
  BOOST_CHECK(rotationallySuperimposable(b, a, Shapes::Shape::Octahedron));

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

  BOOST_CHECK(Temple::testLogicalOperators(c, d));
  BOOST_CHECK(rotationallySuperimposable(c, d, Shapes::Shape::Octahedron));
  BOOST_CHECK(rotationallySuperimposable(d, c, Shapes::Shape::Octahedron));
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
  return Temple::map(
    pairs,
    [](const std::pair<unsigned, unsigned>& v) -> Stereopermutation::Link {
      return std::make_pair(
        Shapes::Vertex(v.first),
        Shapes::Vertex(v.second)
      );
    }
  );
}

void runTestsWithCounts(
  const Shapes::Shape& shape,
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

    // std::cout << stereopermutation.toString() << " unique weights: " << Temple::stringify(unique.weights) << "\n";

    BOOST_CHECK_MESSAGE(
      unique.list.size() == testCase.expectedUnique,
      "Mismatch: Expected " << testCase.expectedUnique
        << " stereopermutations for\n" << stereopermutation.toString() << " in " << Shapes::name(shape) << " shape, got "
        << unique.list.size() << " stereopermutations\n"
    );
  }
}

/* Tetrahedron tests */
BOOST_AUTO_TEST_CASE(MonodentateTetrahedral, *boost::unit_test::label("Stereopermutations")) {
  runTestsWithCounts(
    Shapes::Shape::Tetrahedron,
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
BOOST_AUTO_TEST_CASE(MonodentateSquarePlanar, *boost::unit_test::label("Stereopermutations")) {
  runTestsWithCounts(
    Shapes::Shape::Square,
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
BOOST_AUTO_TEST_CASE(MonodentateOctahedron, *boost::unit_test::label("Stereopermutations")) {
  runTestsWithCounts(
    Shapes::Shape::Octahedron,
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

BOOST_AUTO_TEST_CASE(MultidentateOctahedron, *boost::unit_test::label("Stereopermutations")) {
  runTestsWithCounts(
    Shapes::Shape::Octahedron,
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

BOOST_AUTO_TEST_CASE(OctahedralCaseWithoutTransRemoval, *boost::unit_test::label("Stereopermutations")) {
  const auto shape = Shapes::Shape::Octahedron;
  const auto characters = Characters(6, 'A');
  const auto pairs = PairSet({{0, 1}, {2, 3}, {4, 5}});

  const Stereopermutation stereopermutation {characters, makeLinks(pairs)};
  const auto unique = uniques(stereopermutation, shape);
  BOOST_CHECK_EQUAL(unique.list.size(), 4);
}

BOOST_AUTO_TEST_CASE(numUnlinkedStereopermutationsTest, *boost::unit_test::label("Stereopermutations")) {
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

      const unsigned uniquesCount = uniques(initialStereopermutation, shape).list.size();

      // shapes prediction
      const unsigned chemicalSymmetriesCount = Shapes::Properties::numUnlinkedStereopermutations(
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
