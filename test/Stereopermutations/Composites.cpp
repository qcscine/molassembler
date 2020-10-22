/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include <boost/test/unit_test.hpp>

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Stereopermutation/Composites.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Stringify.h"

#include <iostream>

using namespace Scine::Molassembler;
using namespace Stereopermutations;

namespace {

Shapes::Vertex operator "" _v(unsigned long long v) {
  return Shapes::Vertex(v);
}

std::ostream& operator << (std::ostream& os, const Composite& composite) {
  const unsigned P = composite.allPermutations().size();
  os << P << " permutations:\n";
  for(const auto& p : composite) {
    double dihedral = std::get<2>(p.dihedrals.front());
    os << "- ";
    if(p.alignment == Composite::Alignment::Eclipsed) {
      os << "E ";
    } else if(p.alignment == Composite::Alignment::Staggered) {
      os << "S ";
    }
    os << p.alignedVertices.first << "-" << p.alignedVertices.second << " = " << Temple::Math::toDegrees(dihedral);

    if(p.rankingEquivalentTo) {
      os << "(equivalent to " << p.rankingEquivalentTo->first << "-" << p.rankingEquivalentTo->second << ")";
    }

    os << ", dihedrals " << Temple::stringify(
      Temple::map(p.dihedrals,
        [](const auto& tup) {
          return static_cast<int>(Temple::Math::toDegrees(std::get<2>(tup)));
        }
      )
    );

    os << "\n";
  }

  return os;
}

void expectDominantAngles(
  const Composite& composite,
  const std::vector<double>& angleSet
) {
  for(const auto& permutation : composite) {
    const double dihedral = std::get<2>(permutation.dihedrals.front());

    const auto findIter = std::find_if(
      std::begin(angleSet),
      std::end(angleSet),
      [dihedral](const double v) {
        return Temple::Floating::isCloseAbsolute(v, dihedral, M_PI / 720);
      }
    );

    BOOST_CHECK_MESSAGE(findIter != std::end(angleSet),
      "Could not find dominant angle " << dihedral << " in set of expected dominant angles" << Temple::stringify(angleSet)
    );
  }
}

} // namespace

bool testOrientationState(const Composite::OrientationState& a) {
  // Make a copy and modify that
  auto aCopy = a;

  // Transform and revert the OrientationState
  auto reversionMapping = aCopy.transformToCanonical();
  aCopy.revert(reversionMapping);

  return aCopy == a;
}

// Ensure that transformation and reversion work the way they should
BOOST_AUTO_TEST_CASE(OrientationStateTests, *boost::unit_test::label("Stereopermutations")) {
  for(const auto& shape : Shapes::allShapes) {
    const unsigned S = Shapes::size(shape);
    if(S > 8) {
      continue;
    }

    std::vector<char> maximumAsymmetricCase (S);
    for(unsigned i = 0; i < S; ++i) {
      maximumAsymmetricCase.at(i) = 'A' + i;
    }

    for(Shapes::Vertex i {0}; i < S; ++i ) {
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

BOOST_AUTO_TEST_CASE(CompositePermutationCounts, *boost::unit_test::label("Stereopermutations")) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  struct Args {
    std::pair<Shapes::Shape, Shapes::Shape> shapes;
    std::pair<std::string, std::string> characters;
    std::pair<Shapes::Vertex, Shapes::Vertex> vertices = {0_v, 0_v};
    Composite::Alignment alignment = Composite::Alignment::Eclipsed;

    operator Composite() const {
      auto chars = Temple::mapHomogeneousPairlike(characters,
        [](const auto& str) -> std::vector<char> {
          return {str.begin(), str.end()};
        }
      );

      return Composite {
        Composite::OrientationState {
          shapes.first,
          vertices.first,
          chars.first,
          leftIdentifier
        },
        Composite::OrientationState {
          shapes.second,
          vertices.second,
          chars.second,
          rightIdentifier
        },
        alignment
      };
    }
  };

  struct Expectation {
    unsigned rotations;
    unsigned permutations;
    unsigned axisOrder;
  };

  const std::vector<std::pair<Args, Expectation>> tests {
    {
      Args {
        {Shapes::Shape::Seesaw, Shapes::Shape::Tetrahedron},
        {"ABCD", "ABCD"},
      },
      Expectation {3, 3, 1}
    },
    {
      Args {
        {Shapes::Shape::Octahedron, Shapes::Shape::Octahedron},
        {"ABCDEF", "ABCDEF"},
      },
      Expectation {4, 4, 1}
    },
    {
      Args {
        {Shapes::Shape::Bent, Shapes::Shape::EquilateralTriangle},
        {"AB", "ABC"},
      },
      Expectation {2, 2, 1}
    },
    {
      Args {
        {Shapes::Shape::EquilateralTriangle, Shapes::Shape::EquilateralTriangle},
        {"ABC", "ABC"},
      },
      Expectation {2, 2, 1}
    },
    {
      Args {
        {Shapes::Shape::Bent, Shapes::Shape::Bent},
        {"AB", "AB"},
      },
      Expectation {2, 2, 1}
    },
    {
      Args {
        {Shapes::Shape::Line, Shapes::Shape::Line},
        {"AB", "AB"},
      },
      Expectation {0, 0, 0}
    },
    {
      Args {
        {Shapes::Shape::Octahedron, Shapes::Shape::Octahedron},
        {"ABABCC", "AAAAAA"},
      },
      Expectation {4, 1, 4}
    },
    {
      Args {
        {Shapes::Shape::Octahedron, Shapes::Shape::Octahedron},
        {"ABABCC", "ABABCC"},
      },
      Expectation {4, 2, 2}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::Tetrahedron},
        {"AAAA", "ABCD"},
      },
      Expectation {3, 1, 3}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::VacantTetrahedron},
        {"ABCD", "AAA"},
      },
      Expectation {3, 3, 1}
    },
    {
      Args {
        {Shapes::Shape::EquilateralTriangle, Shapes::Shape::EquilateralTriangle},
        {"ABC", "AAA"},
      },
      Expectation {2, 1, 2}
    },
    {
      Args {
        {Shapes::Shape::EquilateralTriangle, Shapes::Shape::EquilateralTriangle},
        {"AAA", "AAA"},
      },
      Expectation {2, 1, 2}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::EquilateralTriangle},
        {"ABCD", "ABC"},
      },
      Expectation {6, 6, 1}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::EquilateralTriangle},
        {"AAAA", "ABC"},
      },
      Expectation {6, 2, 3}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::EquilateralTriangle},
        {"ABCD", "AAA"},
      },
      Expectation {6, 3, 2}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::EquilateralTriangle},
        {"AAAA", "AAA"},
      },
      Expectation {6, 1, 6}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::EquilateralTriangle},
        {"AAAB", "AAA"},
      },
      Expectation {6, 3, 2}
    },
    {
      Args {
        {Shapes::Shape::EquilateralTriangle, Shapes::Shape::EquilateralTriangle},
        {"BAB", "CAB"},
        {1_v, 2_v}
      },
      Expectation {2, 1, 2}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::Tetrahedron},
        {"ABCD", "AAAA"},
        {0_v, 0_v},
        Composite::Alignment::Staggered
      },
      Expectation {3, 1, 3}
    },
    {
      Args {
        {Shapes::Shape::EquilateralTriangle, Shapes::Shape::EquilateralTriangle},
        {"AAA", "ABC"},
        {0_v, 0_v},
        Composite::Alignment::EclipsedAndStaggered
      },
      Expectation {4, 2, 2}
    },
    {
      Args {
        {Shapes::Shape::EquilateralTriangle, Shapes::Shape::EquilateralTriangle},
        {"ABC", "ABC"},
        {0_v, 0_v},
        Composite::Alignment::EclipsedAndStaggered
      },
      Expectation {4, 4, 1}
    },
    {
      Args {
        {Shapes::Shape::Tetrahedron, Shapes::Shape::Tetrahedron},
        {"AAAA", "AAAB"},
        {0_v, 0_v},
        Composite::Alignment::EclipsedAndStaggered
      },
      Expectation {6, 2, 3}
    },
  };

  for(const auto& pair : tests) {
    const Args& args = pair.first;
    const Expectation& expectation = pair.second;

    // Generate the composite
    const Composite composite = args;

    std::stringstream context;
    context << Shapes::name(args.shapes.first) << "-" << args.characters.first
      << " & " << Shapes::name(args.shapes.second) << "-" << args.characters.second << "\n";
    context << "Composite: " << composite;

    BOOST_TEST_CONTEXT(context.str()) {
      BOOST_CHECK_EQUAL(composite.allPermutations().size(), expectation.rotations);
      BOOST_CHECK_EQUAL(composite.countNonEquivalentPermutations(), expectation.permutations);
      BOOST_CHECK_EQUAL(composite.rotationalAxisSymmetryOrder(), expectation.axisOrder);
    }
  }
}

BOOST_AUTO_TEST_CASE(CompositeAlignment, *boost::unit_test::label("Stereopermutations")) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  Composite a {
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {'A', 'B', 'C', 'D'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {'A', 'B', 'C', 'D'},
      rightIdentifier
    },
    Composite::Alignment::Staggered
  };

  BOOST_CHECK_EQUAL(a.allPermutations().size(), 3);

  Composite b {
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::SquarePyramid,
      0_v,
      {'A', 'B', 'C', 'D', 'E'},
      rightIdentifier
    },
    Composite::Alignment::Staggered
  };

  BOOST_CHECK_EQUAL(b.allPermutations().size(), 3);

  Composite c {
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      4_v,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      2_v,
      {'A', 'B', 'C', 'D', 'E', 'F'},
      rightIdentifier
    }
  };

  BOOST_CHECK_EQUAL(c.allPermutations().size(), 4);
}

/* Expected property of Composite: Rotational variations of ranking characters
 * yield same stereopermutation sequence
 */
BOOST_AUTO_TEST_CASE(CompositeHomomorphisms, *boost::unit_test::label("Stereopermutations")) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  auto testHomomorphism = [](const Composite& a, const Composite& b) {
    // Transform stereopermutations into ranking space
    auto stereopermutationRankingCharacters = [](const Composite& composite) {
      return Temple::map(composite,
        [&](const auto& permutation) {
          return Temple::map(
            permutation.dihedrals,
            [&](Shapes::Vertex i, Shapes::Vertex j, double dihedral) -> std::tuple<char, char, double> {
              return std::make_tuple(
                composite.orientations().first.characters.at(i),
                composite.orientations().second.characters.at(j),
                dihedral
              );
            }
          );
        }
      );
    };

    auto printRankingStereopermutations = [&](const Composite& composite) {
      auto rankingStereo = stereopermutationRankingCharacters(composite);
      for(unsigned i = 0; i < rankingStereo.size(); ++i) {
        std::cout << "Permutation " << i << "\n";
        for(const auto& tup : rankingStereo.at(i)) {
          std::cout << "- " << std::get<0>(tup) << ", " << std::get<1>(tup) << ", " << std::get<2>(tup) << "\n";
        }
      }
      std::cout << "\n";
    };

    const bool pass = Temple::all_of(
      Temple::Adaptors::zip(
        stereopermutationRankingCharacters(a),
        stereopermutationRankingCharacters(b)
      ),
      std::equal_to<>()
    );

    BOOST_CHECK_MESSAGE(pass, "Not all stereopermutation ranking characters are equal!");
    if(!pass) {
      printRankingStereopermutations(a);
      printRankingStereopermutations(b);
    }
  };

  Composite base {
    Composite::OrientationState {
      Shapes::Shape::VacantTetrahedron,
      2_v,
      {{'A', 'A', 'B'}},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      3_v,
      {{'A', 'A', 'B', 'C'}},
      rightIdentifier
    }
  };

  Composite singleRotation {
    Composite::OrientationState {
      Shapes::Shape::VacantTetrahedron,
      0_v,
      {{'B', 'A', 'A'}},
      leftIdentifier
    },
    base.orientations().second
  };
  assert(singleRotation.orientations().second == base.orientations().second);

  testHomomorphism(base, singleRotation);

  Composite multipleRotation {
    singleRotation.orientations().first,
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      1_v,
      {{'B', 'C', 'A', 'A'}},
      rightIdentifier
    }
  };
  assert(multipleRotation.orientations().first == singleRotation.orientations().first);

  testHomomorphism(base, multipleRotation);
}

BOOST_AUTO_TEST_CASE(CompositeCombinedAlignments, *boost::unit_test::label("Stereopermutations")) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  const Composite simplestEclipsed {
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0_v,
      {{'A', 'B'}},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0_v,
      {{'A', 'B'}},
      rightIdentifier
    },
    Composite::Alignment::Eclipsed
  };
  BOOST_CHECK_EQUAL(simplestEclipsed.allPermutations().size(), 2);
  expectDominantAngles(simplestEclipsed, std::vector<double> {{0, M_PI}});

  const Composite simplestStaggered {
    simplestEclipsed.orientations().first,
    simplestEclipsed.orientations().second,
    Composite::Alignment::Staggered
  };
  BOOST_CHECK_EQUAL(simplestStaggered.allPermutations().size(), 2);
  expectDominantAngles(
    simplestStaggered,
    std::vector<double> {{-M_PI / 2, M_PI / 2}}
  );

  const Composite simplestCombined {
    simplestEclipsed.orientations().first,
    simplestEclipsed.orientations().second,
    Composite::Alignment::EclipsedAndStaggered
  };

  BOOST_CHECK_EQUAL(simplestCombined.allPermutations().size(), 4);
  expectDominantAngles(
    simplestCombined,
    std::vector<double> {{-M_PI / 2, 0, M_PI / 2, M_PI}}
  );

  const Composite simplestOffset {
    simplestEclipsed.orientations().first,
    simplestEclipsed.orientations().second,
    Composite::Alignment::BetweenEclipsedAndStaggered
  };
  BOOST_CHECK_EQUAL(simplestOffset.allPermutations().size(), 4);
  expectDominantAngles(
    simplestOffset,
    std::vector<double> {{-3 * M_PI / 4, -M_PI / 4, M_PI / 4, 3 * M_PI / 4}}
  );

  const Composite bothTetrahedral {
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {{'A', 'B', 'C', 'D'}},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {{'A', 'B', 'C', 'D'}},
      rightIdentifier
    },
    Composite::Alignment::EclipsedAndStaggered
  };
  BOOST_CHECK_EQUAL(bothTetrahedral.allPermutations().size(), 6);

  const Composite offsetTetrahedral {
    bothTetrahedral.orientations().first,
    bothTetrahedral.orientations().second,
    Composite::Alignment::BetweenEclipsedAndStaggered
  };
  BOOST_CHECK_EQUAL(offsetTetrahedral.allPermutations().size(), 6);
  expectDominantAngles(
    offsetTetrahedral,
    std::vector<double> {{-5 * M_PI / 6, -3 * M_PI / 6, -M_PI / 6, M_PI / 6, 3 * M_PI / 6, 5 * M_PI /6}}
  );

  /* Eclipsed and staggered stereopermutations for triangle and tetrahedron are
   * identical, so deduplication yields only six
   */
  const Composite bothTriangleTetrahedron {
    Composite::OrientationState {
      Shapes::Shape::EquilateralTriangle,
      0_v,
      {{'A', 'B', 'C'}},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {{'A', 'B', 'C', 'D'}},
      rightIdentifier
    },
    Composite::Alignment::EclipsedAndStaggered
  };
  BOOST_CHECK_EQUAL(bothTriangleTetrahedron.allPermutations().size(), 6);
}
