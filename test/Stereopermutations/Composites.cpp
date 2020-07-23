/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include <boost/test/unit_test.hpp>

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Stereopermutation/Composites.h"
#include "Molassembler/Temple/Adaptors/Zip.h"
#include "Molassembler/Temple/Functional.h"

using namespace Scine::Molassembler;
using namespace Stereopermutations;

namespace {

Shapes::Vertex operator "" _v(unsigned long long v) {
  return Shapes::Vertex(v);
}

} // namespace

bool testOrientationState(Composite::OrientationState a) {
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

BOOST_AUTO_TEST_CASE(CompositeExamples, *boost::unit_test::label("Stereopermutations")) {
  constexpr unsigned leftIdentifier = 0;
  constexpr unsigned rightIdentifier = 1;

  Composite a {
    Composite::OrientationState {
      Shapes::Shape::Seesaw,
      0_v,
      {'A', 'B', 'C', 'D'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {'A', 'B', 'C', 'D'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    a.permutations() == 3u,
    "Expected 3 permutations, got " << a.permutations()
  );

  BOOST_CHECK_MESSAGE(
    !a.isIsotropic(),
    "Expected non-isotropicity for asym seesaw/tetrahedron"
  );

  Composite b {
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

  BOOST_CHECK_MESSAGE(
    b.permutations() == 4u,
    "Expected 4 permutations, got " << b.permutations()
  );

  BOOST_CHECK_MESSAGE(
    !b.isIsotropic(),
    "Expected non-isotropicity for asym oct-oct"
  );

  Composite c {
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::EquilateralTriangle,
      0_v,
      {'A', 'B', 'C'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    c.permutations() == 2,
    "Expected 2 permutations, got " << c.permutations()
  );

  BOOST_CHECK_MESSAGE(
    !c.isIsotropic(),
    "Expected non-isotropicity for asym bent / equilat. triangle"
  );

  Composite d {
    Composite::OrientationState {
      Shapes::Shape::EquilateralTriangle,
      0_v,
      {'A', 'B', 'C'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::EquilateralTriangle,
      0_v,
      {'A', 'B', 'C'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    d.permutations() == 2,
    "Expected 2 permutations, got " << d.permutations()
  );

  BOOST_CHECK_MESSAGE(
    !d.isIsotropic(),
    "Expected non-isotropicity for equilat. triangle/equilat. triangle"
  );

  Composite e {
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Bent,
      0_v,
      {'A', 'B'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    e.permutations() == 2,
    "Expected 2 permutations, got " << e.permutations()
  );

  BOOST_CHECK_MESSAGE(
    !e.isIsotropic(),
    "Expected non-isotropicity for bent/bent"
  );

  Composite f {
    Composite::OrientationState {
      Shapes::Shape::Line,
      0_v,
      {'A', 'B'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Line,
      0_v,
      {'A', 'B'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    f.permutations() == 0,
    "Expected 0 permutations, got " << f.permutations()
  );

  Composite g {
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      0_v,
      {'A', 'B', 'A', 'B', 'C', 'C'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Octahedron,
      0_v,
      {'A', 'A', 'A', 'A', 'A', 'A'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    g.isIsotropic(),
    "Expected isotropicity for matched ranking pairs along bond in oct side"
  );

  Composite h {
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {'A', 'A', 'A', 'A'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {'A', 'B', 'C', 'D'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    h.isIsotropic(),
    "Expected isotropicity for sym. tetr. / asym. tetr. pair"
  );

  Composite i {
    Composite::OrientationState {
      Shapes::Shape::Tetrahedron,
      0_v,
      {'A', 'B', 'C', 'D'},
      leftIdentifier
    },
    Composite::OrientationState {
      Shapes::Shape::VacantTetrahedron,
      0_v,
      {'A', 'A', 'A'},
      rightIdentifier
    }
  };

  BOOST_CHECK_MESSAGE(
    !i.isIsotropic(),
    "Expected non-isotropicity for asym. tetr./sym. vac. tetr pair"
  );
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

  BOOST_CHECK_MESSAGE(
    a.permutations() == 3,
    "Expected three permutations, got " << a.permutations()
  );

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

  BOOST_CHECK_MESSAGE(
    b.permutations() == 3,
    "Expected three permutations, got " << b.permutations()
  );

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

  BOOST_CHECK_MESSAGE(
    c.permutations() == 4,
    "Expected four permutations, got " << c.permutations()
  );
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
        Temple::Functor::map(
          [&](Shapes::Vertex i, Shapes::Vertex j, double dihedral) -> std::tuple<char, char, double> {
            return std::make_tuple(
              composite.orientations().first.characters.at(i),
              composite.orientations().second.characters.at(j),
              dihedral
            );
          }
        )
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
