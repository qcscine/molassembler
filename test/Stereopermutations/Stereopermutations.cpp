/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
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

Rank operator "" _r(unsigned long long i) {
  return Rank(i);
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

    Stereopermutation testStereopermutation {
      Stereopermutation::Occupation {std::vector<unsigned>(S, 0)}
    };
  }
}

inline auto makeOccupation(const std::string& str) {
  return Stereopermutation::occupationFromChars(str);
}

BOOST_AUTO_TEST_CASE(StereopermutationRotation, *boost::unit_test::label("Stereopermutations")) {
  // Octahedron ship-screw like cis-cis-cis
  const Stereopermutation asymm {
    makeOccupation("ABCDEF"),
    {{0_v, 1_v}, {2_v, 4_v}, {3_v, 5_v}}
  };

  std::vector<Shapes::Vertex> C4z {{3_v, 0_v, 1_v, 2_v, 4_v, 5_v}};

  const Stereopermutation expected {
    makeOccupation("DABCEF"),
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
    makeOccupation("AACDBB"),
    {{0_v, 5_v}, {1_v, 4_v}}
  };

  auto isAorB = [](const Rank& rank) -> bool {
    return (rank == 0 || rank == 1);
  };

  const auto allRotations = generateAllRotations(testCase, Shapes::Shape::Octahedron);
  // Make sure the links still refer to A and B 'characters'
  for(const auto& rotation : allRotations) {
    BOOST_CHECK(
      Temple::all_of(
        rotation.links,
        [&](const auto& link) -> bool {
          return (
            isAorB(rotation.occupation.at(link.first))
            && isAorB(rotation.occupation.at(link.second))
          );
        }
      )
    );
  }
}

BOOST_AUTO_TEST_CASE(OctahedralSymmetryCorrectness, *boost::unit_test::label("Stereopermutations")) {
  Stereopermutation octahedralInstance(
    makeOccupation("ABCDEF")
  );

  BOOST_CHECK_EQUAL(
    generateAllRotations(octahedralInstance, Shapes::Shape::Octahedron).size(),
    24 // 4 C4 cases on each of 6 A position selections
  );
}

BOOST_AUTO_TEST_CASE(BugfixTests, *boost::unit_test::label("Stereopermutations")) {
  Stereopermutation a {
    makeOccupation("AAABBB"),
    {{2_v, 3_v}, {1_v, 4_v}, {0_v, 5_v}}
  };
  Stereopermutation b {
    makeOccupation("AABABB"),
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
    makeOccupation("AAABBB"),
    {{2_v, 3_v}, {1_v, 4_v}, {0_v, 5_v}}
  };
  Stereopermutation d {
    makeOccupation("AAABBB"),
    {{1_v, 4_v}, {2_v, 3_v}, {0_v, 5_v}}
  };

  BOOST_CHECK(Temple::testLogicalOperators(c, d));
  BOOST_CHECK(rotationallySuperimposable(c, d, Shapes::Shape::Octahedron));
  BOOST_CHECK(rotationallySuperimposable(d, c, Shapes::Shape::Octahedron));
}

using Occupation = Stereopermutation::Occupation;
using PairSet = std::vector<
  std::pair<unsigned, unsigned>
>;

struct TestCase {
  Occupation occupation;
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
      testCase.occupation,
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
      {makeOccupation("AAAA"), {}, 1},
      {makeOccupation("AAAB"), {}, 1},
      {makeOccupation("AABB"), {}, 1},
      {makeOccupation("AABC"), {}, 1},
      {makeOccupation("ABCD"), {}, 2},
    }
  );
}

/* Square Planar tests */
BOOST_AUTO_TEST_CASE(MonodentateSquarePlanar, *boost::unit_test::label("Stereopermutations")) {
  runTestsWithCounts(
    Shapes::Shape::Square,
    {
      {makeOccupation("AAAA"), {}, 1},
      {makeOccupation("AAAB"), {}, 1},
      {makeOccupation("AABB"), {}, 2},
      {makeOccupation("AABC"), {}, 2},
      {makeOccupation("ABCD"), {}, 3},
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
      {makeOccupation("AAAAAA"), {}, 1},
      {makeOccupation("AAAAAB"), {}, 1},
      {makeOccupation("AAAABB"), {}, 2},
      {makeOccupation("AAABBB"), {}, 2},
      {makeOccupation("AAAABC"), {}, 2},
      {makeOccupation("AAABBC"), {}, 3},
      {makeOccupation("AABBCC"), {}, 6},
      {makeOccupation("AAABCD"), {}, 5},
      {makeOccupation("AABBCD"), {}, 8},
      {makeOccupation("AABCDE"), {}, 15},
      {makeOccupation("ABCDEF"), {}, 30},
    }
  );
}

BOOST_AUTO_TEST_CASE(MultidentateOctahedron, *boost::unit_test::label("Stereopermutations")) {
  runTestsWithCounts(
    Shapes::Shape::Octahedron,
    {
      { // M(A-A)_3
        makeOccupation("AAAAAA"),
        PairSet({{0, 1}, {2, 3}, {4, 5}}),
        2
        /* NOTE: besides the two cis-cis-cis enantiomers, there are
         * cis-cis-trans and trans-trans-trans isomers, these are removed by
         * default!
         */
      },
      { // M(A-B)_3
        makeOccupation("ABABAB"),
        PairSet({{0, 1}, {2, 3}, {4, 5}}),
        4
      },
      { // M(A-B)_2 CD
        makeOccupation("ABABCD"),
        PairSet({{0, 1}, {2, 3}}),
        11
      },
      { // M(A-A)(B-C)DE
        makeOccupation("AABCDE"),
        PairSet({{0, 1}, {2, 3}}),
        10
      },
      { // M(A-B)(C-D)EF
        makeOccupation("ABCDEF"),
        PairSet({{0, 1}, {2, 3}}),
        20
      },
      { // M(A-B-A)CDE
        makeOccupation("ABACDE"),
        PairSet({{0, 1}, {1, 2}}),
        9
      },
      { // M(A-B-C)_2
        makeOccupation("ABCABC"),
        PairSet({{0, 1}, {1, 2}, {3, 4}, {4, 5}}),
        11
      },
      { // M(A-B-B-A)CD
        makeOccupation("ABBACD"),
        PairSet({{0, 1}, {1, 2}, {2, 3}}),
        7
      },
      { // M(A-B-C-B-A)D
        makeOccupation("ABCBAD"),
        PairSet({{0, 1}, {1, 2}, {2, 3}, {3, 4}}),
        7
      }
    }
  );
}

BOOST_AUTO_TEST_CASE(OctahedralCaseWithoutTransRemoval, *boost::unit_test::label("Stereopermutations")) {
  const auto shape = Shapes::Shape::Octahedron;
  const auto occupation = makeOccupation("AAAAAA");
  const auto pairs = PairSet({{0, 1}, {2, 3}, {4, 5}});

  const Stereopermutation stereopermutation {occupation, makeLinks(pairs)};
  const auto unique = uniques(stereopermutation, shape);
  BOOST_CHECK_EQUAL(unique.list.size(), 4);
}

BOOST_AUTO_TEST_CASE(numUnlinkedStereopermutationsTest, *boost::unit_test::label("Stereopermutations")) {
#ifdef NDEBUG
  constexpr unsigned sizeLimit = 8;
#else
  constexpr unsigned sizeLimit = 6;
#endif


  // Crosscheck number of unlinked stereopermutations with shapes
  for(const Shapes::Shape shape : Shapes::allShapes) {
    const unsigned S = Shapes::size(shape);
    if(S > sizeLimit) {
      continue;
    }

    for(unsigned nIdentical = 1; nIdentical < S; ++nIdentical) {
      // Populate characters for construction of a Stereopermutation
      std::string characters;
      for(unsigned i = 0; i < S; ++i) {
        characters += static_cast<char>('A' + std::max(
          0,
          static_cast<int>(i) - static_cast<int>(nIdentical) + 1
        ));
      }

      const unsigned count = std::count(
        std::begin(characters),
        std::end(characters),
        'A'
      );

      BOOST_REQUIRE_EQUAL(count, nIdentical);

      const Stereopermutation initialStereopermutation {makeOccupation(characters)};

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
