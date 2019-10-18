/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Central symmetry data class definitions
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_PRIMITIVES_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_PRIMITIVES_H

#include "temple/constexpr/Vector.h"
#include "temple/constexpr/TupleType.h"

#include "chemical_symmetries/Shapes.h"
#include "chemical_symmetries/PointGroups.h"
#include "chemical_symmetries/CompileTimeOptions.h"
#include "chemical_symmetries/AngleLookup.h"

namespace Scine {

namespace Symmetry {

namespace concepts {

/**
 * @brief Concept checking class for symmetry classes
 *
 * @tparam T Class to check the concept for
 *
 * All members must be `static constexpr` (`const` is then implicit)
 *
 * @par Member requirements (Type and name)
 * - `Shape shape`: Enum member specifying the symmetry name
 * - `unsigned size`: Number of symmetry positions in the symmetry
 * - `char* stringName`: Human readable string of the name
 * - `double(const unsigned, const unsigned) angleFunction`: Angle in radians
 *   between symmetry position indices
 * - `std::array<temple::Vector, size> coordinate`: Origin-centered normalized
 *   position vectors of the symmetry positions
 * - `std::array< std::array<unsigned, size>, ?> rotations`: Spatial rotations
 *   represented as index permutations between symmetry positions. A minimal
 *   set that combined can generate all rotations.
 * - `std::array< std::array<unsigned, 4>, ?> tetrahedra`: A list of tetrahedra
 *   definitions whose signed volumes can uniquely identify a maximally
 *   asymmetric set of ligands
 * - `std::array<unsigned, 0 or size> mirror`: A mirroring symmetry element or
 *   an empty array if the symmetry cannot be chiral.
 */
template<typename T>
struct ShapeClass : std::integral_constant<bool,
  (
    std::is_same<Shape, std::decay_t<decltype(T::shape)>>::value
    && std::is_same<PointGroup, std::decay_t<decltype(T::pointGroup)>>::value
    && std::is_same<unsigned, std::decay_t<decltype(T::size)>>::value
    && std::is_same<const char*, std::decay_t<decltype(T::stringName)>>::value
    && std::is_same<double, decltype(T::angleFunction(0u, 1u))>::value
    && std::is_same<
      temple::Vector,
      temple::getValueType<decltype(T::coordinates)>
    >::value
    && T::coordinates.size() == T::size
    && std::is_same<
      std::array<unsigned, T::size>,
      temple::getValueType<decltype(T::rotations)>
    >::value
    && std::is_same<
      std::array<unsigned, 4>,
      temple::getValueType<decltype(T::tetrahedra)>
    >::value
    && std::is_same<unsigned, temple::getValueType<decltype(T::mirror)>>::value
    && (T::mirror.size() == 0 || T::mirror.size() == T::size)
  )
> {};

} // namespace concepts

//! A placeholder value for constexpr tetrahedra specification of origin
constexpr unsigned ORIGIN_PLACEHOLDER = std::numeric_limits<unsigned>::max();

/*!
 * @brief All symmetry data classes and some minor helper functions
 *
 * Each symmetry data class must fulfill concepts::ShapeClass
 */
namespace data {

/*!
 * @brief Line shape
 *
 * @verbatim
 *
 *  0 – (_) – 1
 *
 * @endverbatim
 */
struct Line {
  static constexpr Shape shape = Shape::Line;
  static constexpr PointGroup pointGroup = PointGroup::Cinfv;
  static constexpr unsigned size = 2;
  static constexpr char stringName[] = "line";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return M_PI;
  }
  static constexpr std::array<temple::Vector, 2> coordinates {{
    { 1 , 0, 0 },
    { -1, 0, 0 }
  }};
  static constexpr std::array<
    std::array<unsigned, 2>,
    1
  > rotations {{
    {{1, 0}}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
  static constexpr std::array<unsigned, 0> mirror {};
};

/*!
 * @brief Bent symmetry
 *
 * @verbatim
 *
 *  1
 *   \
 *    (_) – 0
 *
 * @endverbatim
 */
struct Bent {
  static constexpr Shape shape = Shape::Bent;
  static constexpr PointGroup pointGroup = PointGroup::C2v;
  static constexpr unsigned size = 2;
  static constexpr char stringName[] = "bent";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    /* subject to a lot of variation, between 90 and 109 degrees pursuant to
     * english wikipedia, using experimental data here to improve instances
     * of this geometry on e.g. O center would a big improvement to DG runs
     */
    if(a == b) {
      return 0;
    }

    return temple::Math::toRadians<double>(107);
  }
  static constexpr std::array<temple::Vector, 2> coordinates {{
    {1., 0., 0.},
    {-0.292372, 0.956305, 0.}
  }};
  static constexpr std::array<
    std::array<unsigned, 2>,
    1
  > rotations {{
    {{1, 0}}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
  static constexpr std::array<unsigned, 0> mirror {};
};

/*!
 * @brief Trigonal planar symmetry
 *
 * @verbatim
 *
 *     0
 *     |
 *    (_)
 *   /   \
 *  1     2
 *
 * @endverbatim
 *
 * This character art is not quite ideal since the angles are thoroughly
 * misrepresented, but all positions including the central atom are in one
 * plane. The angles are idealized as 120°.
 */
struct EquilateralTriangle {
  static constexpr Shape shape = Shape::EquilateralTriangle;
  static constexpr PointGroup pointGroup = PointGroup::D3h;
  static constexpr unsigned size = 3;
  static constexpr char stringName[] = "triangle";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return temple::Math::toRadians<double>(120);
  }
  static constexpr std::array<temple::Vector, 3> coordinates {{
      {1, 0, 0},
      {-0.5, 0.866025, 0},
      {-0.5, -0.866025, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    2
  > rotations {{
    {{1, 2, 0}}, // C3
    {{0, 2, 1}} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
  static constexpr std::array<unsigned, 0> mirror {};
};

/*!
 * @brief A Tetrahedron symmetry missing a ligand
 *
 * This symmetry is widely called trigonal pyramidal, but this is clearly a
 * misnomer. Trigonal pyramidal should denote the symmetry that is trigonal
 * planar plus an axial ligand, one short of trigonal bipyramidal.
 *
 * @verbatim
 *
 *     (_)
 *    /  \ °2
 *   0    1
 *
 * Where /, \ denote in front of plane bonds, ° a behind the plane bond.
 *
 * @endverbatim
 */
struct ApicalTrigonalPyramid {
  static constexpr Shape shape = Shape::ApicalTrigonalPyramid;
  static constexpr PointGroup pointGroup = PointGroup::C3v;
  static constexpr unsigned size = 3;
  static constexpr char stringName[] = "apical trigonal pyramid";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return temple::Math::toRadians<double>(107.5);
  }
  static constexpr std::array<temple::Vector, 3> coordinates {{
    {0, -0.366501, 0.930418},
    {0.805765, -0.366501, -0.465209},
    {-0.805765, -0.366501, -0.465209}
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    1
  > rotations {{
    {{2, 0, 1}} // C3
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{ORIGIN_PLACEHOLDER, 0, 1, 2}}
  }};
  static constexpr std::array<unsigned, 3> mirror {{0, 2, 1}};
};

/*!
 * @brief A T-shaped symmetry
 *
 * @verbatim
 *
 * 0 – (_) – 2
 *      |
 *      1
 *
 * @endverbatim
 */
struct T {
  static constexpr Shape shape = Shape::T;
  static constexpr PointGroup pointGroup = PointGroup::C2v;
  static constexpr unsigned size = 3;
  static constexpr char stringName[] = "T-shaped";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if((a + b) % 2 == 1) {
      return M_PI / 2;
    }

    return M_PI;
  }
  static constexpr std::array<temple::Vector, 3> coordinates {{
    {-1, -0, -0},
    {0, 1, 0},
    {1, 0, 0},
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    1
  > rotations {{
    {{2, 1, 0}} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
  static constexpr std::array<unsigned, 0> mirror {};
};

/*!
 * @brief A tetrahedral symmetry
 *
 * @verbatim
 *
 *      0
 *      |
 *     (_)
 *    /  \ °3
 *   1    2
 *
 * Where /, \ denote in front of plane bonds, ° a behind the plane bond.
 *
 * @endverbatim
 */
struct Tetrahedron {
  static constexpr Shape shape = Shape::Tetrahedron;
  static constexpr PointGroup pointGroup = PointGroup::Td;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "tetrahedron";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return temple::Math::toRadians<double>(109.5);
  }
  static constexpr std::array<temple::Vector, 4> coordinates {{
    {0, 1, 0},
    {0, -0.333807, 0.942641},
    {0.816351, -0.333807, -0.471321},
    {-0.816351, -0.333807, -0.471321}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > rotations {{
    {{0, 3, 1, 2}}, // C4, 1
    {{2, 1, 3, 0}}, // C4, 2
    {{3, 0, 2, 1}}, // C4, 3
    {{1, 2, 0, 3}}  // C4, 4
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{0, 1, 2, 3}}
  }};
  static constexpr std::array<unsigned, 4> mirror {{0, 2, 1, 3}};
};

/*!
 * @brief A square planar symmetry
 *
 * @verbatim
 *
 *   3   2
 *    \_/
 *    (_) <- central atom
 *    / \
 *   0   1
 *
 * @endverbatim
 */
struct Square {
  static constexpr Shape shape = Shape::Square;
  static constexpr PointGroup pointGroup = PointGroup::D4h;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "square";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if((a + b) % 2 == 1) {
      // this expression indicates cis
      return M_PI / 2;
    }

    // leftover case is trans
    return M_PI;
  }
  static constexpr std::array<temple::Vector, 4> coordinates {{
    {1, 0, 0},
    {0, 1, 0},
    {-1, -0, -0},
    {-0, -1, -0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > rotations {{
    {{3, 0, 1, 2}}, // C4
    {{1, 0, 3, 2}}, // C2
    {{3, 2, 1, 0}}  // C2'
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
  static constexpr std::array<unsigned, 0> mirror {};
};

/*!
 * @brief A disphenoid edge-centered / seesaw shape
 *
 * @verbatim
 *
 *   0 – (_) – 3
 *       / :
 *      1   2
 *
 * @endverbatim
 */
struct Disphenoid {
  static constexpr Shape shape = Shape::Disphenoid;
  static constexpr PointGroup pointGroup = PointGroup::C2v;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "disphenoid";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    const auto& smaller = std::min(a, b);
    const auto& larger = std::max(a, b);
    if(smaller == 0 && larger == 3) {
      return M_PI;
    }

    if(smaller == 1 && larger == 2) {
      return temple::Math::toRadians<double>(120);
    }

    return M_PI / 2;
  }
  static constexpr std::array<temple::Vector, 4> coordinates {{
    {0, 1, 0},
    {1, 0, 0},
    {-0.5, 0, -0.866025},
    {-0, -1, -0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > rotations {{
    {{3, 2, 1, 0}} // C2
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{0, 1, 2, 3}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, ORIGIN_PLACEHOLDER, 1, 2}},
    {{ORIGIN_PLACEHOLDER, 3, 1, 2}},
  }};
#endif

  static constexpr std::array<unsigned, 4> mirror {{0, 2, 1, 3}};
};

/*!
 * @brief A trigonal pyramidal symmetry = trig. pl. + an axial ligand
 *
 * A trigonal planar symmetry + an axial ligand.
 *
 * @verbatim
 *
 * Viewed from the top of the pyramid. The central atom is ( ), 3 is apical
 *
 *       3
 *       |  2
 *       | :    <- behind view plane
 *   0--(_)
 *         \    <- in front of view plane
 *          1
 *
 * @endverbatim
 */
struct TrigonalPyramid {
  static constexpr Shape shape = Shape::TrigonalPyramid;
  static constexpr PointGroup pointGroup = PointGroup::C3v;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "trigonal pyramid";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if(std::max(a, b) != 3) {
      // -> smaller < 2, this means either 0,1 0,2 1,2 axial
      return temple::Math::toRadians<double>(120);
    }

    // -> smaller < 3, this means {1,2,3}, 3
    return M_PI / 2;
  }
  static constexpr std::array<temple::Vector, 4> coordinates {{
    {1, 0, 0},
    {-0.5, 0.866025, 0},
    {-0.5, -0.866025, 0},
    {0, 0, 1}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > rotations {{
    {{2, 0, 1, 3}} // C3
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{0, 1, 3, 2}}
  }};
  static constexpr std::array<unsigned, 4> mirror {{0, 2, 1, 3}};
};

/*!
 * @brief A square pyramidal symmetry
 *
 * @verbatim
 *
 *      4
 *   3  |  2
 *    : | :    <- behind view plane
 *     (_)
 *    /   \    <- in front of view plane
 *   0     1
 *
 * @endverbatim
 */
struct SquarePyramid {
  static constexpr Shape shape = Shape::SquarePyramid;
  static constexpr PointGroup pointGroup = PointGroup::C4v;
  static constexpr unsigned size = 5;
  static constexpr char stringName[] = "square pyramid";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if(a == 4 || b == 4) { // all bonds to axial ligand are 90°
      return M_PI / 2;
    }

    if((a + b) % 2 == 0) { // 0 + 2 or 1 + 3 are trans
      return M_PI;
    }

    // rest are cis
    return M_PI / 2;
  }
  static constexpr std::array<temple::Vector, 5> coordinates {{
    {1, 0, 0},
    {0, 1, 0},
    {-1, -0, -0},
    {-0, -1, -0},
    {0, 0, 1}
  }};
  static constexpr std::array<
    std::array<unsigned, 5>,
    1
  > rotations {{
    {{3, 0, 1, 2, 4}} // C4
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, 1, 4, 2}},
    {{0, 3, 2, 4}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
    {{0, 1, 4, ORIGIN_PLACEHOLDER}},
    {{1, 2, 4, ORIGIN_PLACEHOLDER}},
    {{2, 3, 4, ORIGIN_PLACEHOLDER}},
    {{3, 0, 4, ORIGIN_PLACEHOLDER}}
  }};
#endif

  static constexpr std::array<unsigned, 5> mirror {{1, 0, 3, 2, 4}};
};

/*!
 * @brief A trigonal bipyramidal symmetry
 *
 * @verbatim
 *
 * Viewed from the top of the pyramid. The central atom is ( ), 3 and 4
 * are axial.
 *
 *       3
 *       |  2
 *       | :    <- behind view plane
 *   0--(_)
 *       | \    <- in front of view plane
 *       |  1
 *       4
 *
 * @endverbatim
 */
struct TrigonalBipyramid {
  static constexpr Shape shape = Shape::TrigonalBipyramid;
  static constexpr PointGroup pointGroup = PointGroup::D3h;
  static constexpr unsigned size = 5;
  static constexpr char stringName[] = "trigonal bipyramid";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    unsigned smaller = std::min(a, b), larger = std::max(a, b);
    if(larger < 3) {
      // -> smaller < 2, this means either 0,1 0,2 1,2 axial
      return 2 * M_PI / 3;
    }

    if(larger == 3) {
      // -> smaller < 3, this means {1,2,3}, 3
      return M_PI / 2;
    }

    if(smaller < 3) {
      // now, larger must be 4 (process of elimination), so if a is not 3:
      return M_PI / 2;
    }

    // only case left: 3,4
    return M_PI;
  }
  static constexpr std::array<temple::Vector, 5> coordinates {{
    {1, 0, 0},
    {-0.5, 0.866025, 0},
    {-0.5, -0.866025, 0},
    {0, 0, 1},
    {-0, -0, -1}
  }};
  static constexpr std::array<
    std::array<unsigned, 5>,
    4
  > rotations {{
    {{2, 0, 1, 3, 4}}, // C3
    {{0, 2, 1, 4, 3}}, // C2 on 0
    {{2, 1, 0, 4, 3}}, // C2 on 1
    {{1, 0, 2, 4, 3}} // C2 on 2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
      {{0, 1, 3, 2}},
      {{0, 1, 2, 4}}
  }};
  static constexpr std::array<unsigned, 5> mirror {{0, 2, 1, 3, 4}};
};

/*!
 * @brief A pentagonal planar shape
 *
 * @verbatim
 *
 * All in plane:
 *
 *      0
 *  1.  |  .4
 *    °(_)°
 *    /   \
 *   2     3
 *
 * @endverbatim
 */
struct Pentagon {
  static constexpr Shape shape = Shape::Pentagon;
  static constexpr PointGroup pointGroup = PointGroup::D5h;
  static constexpr unsigned size = 5;
  static constexpr char stringName[] = "pentagon";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * temple::Math::toRadians<double>(72);
  }
  static constexpr std::array<temple::Vector, 5> coordinates {{
    {1, 0, 0},
    {0.309017, 0.951057, 0},
    {-0.809017, 0.587785, 0},
    {-0.809017, -0.587785, 0},
    {0.309017, -0.951057, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 5>,
    2
  > rotations {{
    {{4, 0, 1, 2, 3}}, // C5
    {{0, 4, 3, 2, 1}} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
  static constexpr std::array<unsigned, 0> mirror {};
};

/*!
 * @brief An octahedral symmetry.
 *
 * @verbatim
 *
 * The central atom is ( ), 4 and 5 are axial, the rest equatorial.
 *
 *       4
 *    3  |  2
 *     : | :
 *      (_)
 *     / | \
 *    0  |  1
 *       5
 *
 * Where /, \ denote bonds in front of the view plane, : denotes bonds
 * behind the view plane.
 *
 * @endverbatim
 */
struct Octahedron {
  static constexpr Shape shape = Shape::Octahedron;
  static constexpr PointGroup pointGroup = PointGroup::Oh;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "octahedron";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if(
      (
        std::max(a, b) < 4 // if the largest is < 4, then equatorial
        && (a + b) % 2 == 0 // this gives trans eq ligands
      ) || std::min(a, b) == 4 // this indicates 4,5 (axial trans)
    ) {
      return M_PI;
    }

    return M_PI / 2;
  }
  static constexpr std::array<temple::Vector, 6> coordinates {{
    {1, 0, 0},
    {0, 1, 0},
    {-1, -0, -0},
    {-0, -1, -0},
    {0, 0, 1},
    {-0, -0, -1}
  }};
  static constexpr std::array<
    std::array<unsigned, 6>,
    3
  > rotations {{
    {{3, 0, 1, 2, 4, 5}}, // vertical C4
    {{0, 5, 2, 4, 1, 3}}, // horizontal C4
    {{4, 1, 5, 3, 2, 0}} // horizontal C4'
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
    {{3, 0, 4, 5}},
    {{0, 1, 4, 5}},
    {{1, 2, 4, 5}},
    {{2, 3, 4, 5}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    8
  > tetrahedra {{
    {{3, 0, 4, ORIGIN_PLACEHOLDER}},
    {{0, 1, 4, ORIGIN_PLACEHOLDER}},
    {{1, 2, 4, ORIGIN_PLACEHOLDER}},
    {{2, 3, 4, ORIGIN_PLACEHOLDER}},
    {{3, 0, ORIGIN_PLACEHOLDER, 5}},
    {{0, 1, ORIGIN_PLACEHOLDER, 5}},
    {{1, 2, ORIGIN_PLACEHOLDER, 5}},
    {{2, 3, ORIGIN_PLACEHOLDER, 5}}
  }};
#endif
  static constexpr std::array<unsigned, 6> mirror {{1, 0, 3, 2, 4, 5}};
};

/*!
 * @brief A trigonal prismatic symmetry
 *
 * @verbatim
 *
 *   3  4  5
 *    : | :
 *     (_)
 *    : | :
 *   0  1  2
 *
 * Where /, \ denote bonds in front of the view plane, : denotes bonds
 * behind the view plane.
 *
 * @endverbatim
 *
 * Angles
 *  0-1, 0-2 -> 86°
 *  0-3 -> 76°
 *  0-4, 0-5 -> 134°
 *
 * From [W(CH3)6], Haaland, Hammel, Rypdal, Volden, J. Am. Chem. Soc. 1990
 */
struct TrigonalPrism {
  static constexpr Shape shape = Shape::TrigonalPrism;
  static constexpr PointGroup pointGroup = PointGroup::D3h;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "trigonal prism";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    // Between plane symmetric
    if(std::min(a - b, b - a) == 3) {
      return temple::Math::toRadians<double>(76);
    }

    // In plane triangle
    if(
      (a < 3 && b < 3)
      || (a >= 3 && b >= 3)
    ) {
      return temple::Math::toRadians<double>(86);
    }

    // Between plane asymmetric
    return temple::Math::toRadians<double>(134);
  }
  static constexpr std::array<temple::Vector, 6> coordinates {{
    {0.788011, 0, -0.615661},
    {-0.394005, 0.682437, -0.615661},
    {-0.394005, -0.682437, -0.615661},
    {0.788011, 0, 0.615661},
    {-0.394005, 0.682437, 0.615661},
    {-0.394005, -0.682437, 0.615661}
  }};
  static constexpr std::array<
    std::array<unsigned, 6>,
    2
  > rotations {{
    {{2, 0, 1, 5, 3, 4}}, // C3 axial
    {{5, 4, 3, 2, 1, 0}} // C2 betw. 1, 4
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{ORIGIN_PLACEHOLDER, 0, 1, 2}},
    {{3, ORIGIN_PLACEHOLDER, 4, 5}}
  }};
  static constexpr std::array<unsigned, 6> mirror {{2, 1, 0, 5, 4, 3}};
};

/*!
 * @brief A pentagonal pyramidal symmetry
 *
 * @verbatim
 *
 *      0
 *  1.  |  .4
 *    °(5)°
 *    :   :
 *   2     3
 *
 * 0-4 are in plane,
 * 5 is above plane,
 * ( ) signifies the central atom beneath it
 *
 * @endverbatim
 */
struct PentagonalPyramid {
  static constexpr Shape shape = Shape::PentagonalPyramid;
  static constexpr PointGroup pointGroup = PointGroup::C5v;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "pentagonal pyramid";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if(a == 5 || b == 5) {
      return M_PI / 2;
    }

    // remainder are identical to Pentagon
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * temple::Math::toRadians<double>(72);
  }
  static constexpr std::array<temple::Vector, 6> coordinates {{
    {1, 0, 0},
    {0.309017, 0.951057, 0},
    {-0.809017, 0.587785, 0},
    {-0.809017, -0.587785, 0},
    {0.309017, -0.951057, 0},
    {0, 0, 1}
  }};
  static constexpr std::array<
    std::array<unsigned, 6>,
    1
  > rotations {{
    {{4, 0, 1, 2, 3, 5}} // C5 axial
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > tetrahedra {{
    {{0, 1, 5, 2}},
    {{2, 3, 5, 4}},
    {{4, 5, ORIGIN_PLACEHOLDER, 0}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    5
  > tetrahedra {{
    {{0, ORIGIN_PLACEHOLDER, 1, 5}},
    {{1, ORIGIN_PLACEHOLDER, 2, 5}},
    {{2, ORIGIN_PLACEHOLDER, 3, 5}},
    {{3, ORIGIN_PLACEHOLDER, 4, 5}},
    {{4, ORIGIN_PLACEHOLDER, 0, 5}}
  }};
#endif
  static constexpr std::array<unsigned, 6> mirror {{0, 4, 3, 2, 1, 5}};
};

/**
 * @brief Hexagonal planar shape
 */
struct Hexagon {
  static constexpr Shape shape = Shape::Hexagon;
  static constexpr PointGroup pointGroup = PointGroup::D6h;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "hexagon";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 6, 6 - absDiff)
    ) * temple::Math::toRadians<double>(60);
  }
  static constexpr std::array<temple::Vector, 6> coordinates {{
    { 1.000000,  0.000000,  0.000000},
    { 0.500000,  0.866025,  0.000000},
    {-0.500000,  0.866025,  0.000000},
    {-1.000000,  0.000000,  0.000000},
    {-0.500000, -0.866025,  0.000000},
    { 0.500000, -0.866025,  0.000000}
  }};
  static constexpr std::array<
    std::array<unsigned, 6>,
    2
  > rotations {{
    {{5, 0, 1, 2, 3, 4}}, // C6
    {{0, 5, 4, 3, 2, 1}} // C2
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
  static constexpr std::array<unsigned, 0> mirror {};
};

/*!
 * @brief A pentagonal bipyramidal symmetry
 *
 * @verbatim
 *
 * 3, 5, (_) and 6 in plane, 1 and 2 in front, 0 and 4 in back
 *
 *      5
 *  0_  | .4
 *   _:(_) – 3
 *  1   |°·2
 *      6
 *
 * 0-4 are equatorial,
 * 5,6 are axial
 *
 * @endverbatim
 */
struct PentagonalBipyramid {
  static constexpr Shape shape = Shape::PentagonalBipyramid;
  static constexpr PointGroup pointGroup = PointGroup::D5h;
  static constexpr unsigned size = 7;
  static constexpr char stringName[] = "pentagonal bipyramid";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if(a + b == 11) {
      return M_PI; // trans 5,6
    }

    if(temple::Math::XOR(a > 4, b > 4)) {
      return M_PI / 2; // any angle to axial index
    }

    // remainder are equatorial angles, like Pentagon
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * temple::Math::toRadians<double>(72);
  }
  static constexpr std::array<temple::Vector, 7> coordinates {{
    {1, 0, 0},
    {0.309017, 0.951057, 0},
    {-0.809017, 0.587785, 0},
    {-0.809017, -0.587785, 0},
    {0.309017, -0.951057, 0},
    {0, 0, 1},
    {-0, -0, -1}
  }};
  static constexpr std::array<
    std::array<unsigned, 7>,
    2
  > rotations {{
    {{4, 0, 1, 2, 3, 5, 6}}, // C5 axial
    {{1, 0, 4, 3, 2, 6, 5}} // C2 equatorial on 3
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    5
  > tetrahedra {{
      {{0, 1, 5, 6}},
      {{1, 2, 5, 6}},
      {{2, 3, 5, 6}},
      {{3, 4, 5, 6}},
      {{4, 0, 5, 6}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    10
  > tetrahedra {{
    {{0, 1, 5, ORIGIN_PLACEHOLDER}},
    {{1, 2, 5, ORIGIN_PLACEHOLDER}},
    {{2, 3, 5, ORIGIN_PLACEHOLDER}},
    {{3, 4, 5, ORIGIN_PLACEHOLDER}},
    {{4, 0, 5, ORIGIN_PLACEHOLDER}},
    {{0, 1, ORIGIN_PLACEHOLDER, 6}},
    {{1, 2, ORIGIN_PLACEHOLDER, 6}},
    {{2, 3, ORIGIN_PLACEHOLDER, 6}},
    {{3, 4, ORIGIN_PLACEHOLDER, 6}},
    {{4, 0, ORIGIN_PLACEHOLDER, 6}}
  }};
#endif
  static constexpr std::array<unsigned, 7> mirror {{0, 4, 3, 2, 1, 5, 6}};
};

/*!
 * @brief A capped octahedron shape
 */
struct CappedOctahedron {
  static constexpr Symmetry::Shape shape = Symmetry::Shape::CappedOctahedron;
  static constexpr unsigned size = 7;
  static constexpr char stringName[] = "capped octahedron";
  /*! [V(CO)7]+ in C3v
   *
   * Jay W. Dicke, Nathan J. Stibrich, Henry F. Schaefer,
   * V(CO)7+: A capped octahedral structure completes the 18-electron rule,
   * Chemical Physics Letters, Volume 456, Issues 1–3, 2008.
   */
  static constexpr std::array<temple::Vector, 7> coordinates {{
    { 0.000000,  0.000000,  1.000000},
    { 0.956305,  0.000000,  0.292372},
    {-0.478152,  0.828184,  0.292372},
    {-0.478152, -0.828184,  0.292372},
    { 0.400888,  0.694358, -0.597625},
    {-0.801776,  0.000000, -0.597625},
    { 0.400888, -0.694358, -0.597625}
  }};
  static constexpr auto angleLookupTable = temple::makeUpperTriangularMatrix(
    detail::makeArray<size>(coordinates)
  );
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return angleLookupTable.at(
      std::min(a, b),
      std::max(a, b)
    );
  }
  static constexpr std::array<
    std::array<unsigned, 7>,
    1
  > rotations {{
    {{0, 3, 1, 2, 6, 4, 5}} // C3 axial
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, 1, 2, 3}},
    {{0, 4, 5, 6}}
  }};
  static constexpr std::array<unsigned, 7> mirror {{0, 3, 2, 1, 6, 5, 4}};
};

/*!
 * @brief A capped trigonal prism shape
 */
struct CappedTrigonalPrism {
  static constexpr Symmetry::Shape shape = Symmetry::Shape::CappedTrigonalPrism;
  static constexpr unsigned size = 7;
  static constexpr char stringName[] = "capped trigonal prism";
  //! [V(CO)7]+ in C2v, same as from CappedOctahedron
  static constexpr std::array<temple::Vector, 7> coordinates {{
    { 0.000000,  0.000000,  1.000000},
    { 0.990268,  0.000000,  0.139173},
    { 0.000000,  0.990268,  0.139173},
    {-0.990268,  0.000000,  0.139173},
    {-0.000000, -0.990268,  0.139173},
    { 0.414628,  0.414628, -0.810042},
    {-0.414628, -0.414628, -0.810042}
  }};
  static constexpr auto angleLookupTable = temple::makeUpperTriangularMatrix(
    detail::makeArray<size>(coordinates)
  );
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return angleLookupTable.at(
      std::min(a, b),
      std::max(a, b)
    );
  }
  static constexpr std::array<
    std::array<unsigned, 7>,
    1
  > rotations {{
    {{0, 3, 4, 1, 2, 6, 5}} // C2 axial
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, 1, 2, 5}},
    {{0, 3, 4, 6}}
  }};
  static constexpr std::array<unsigned, 7> mirror {{0, 2, 1, 4, 3, 5, 6}};
};

/*!
 * @brief A square antiprismatic symmetry
 *
 * @verbatim
 *
 * Two representations, one oblique, the other more helpful.
 * The first is a side-on view. 0 is mostly hidden by 1, 3 mostly hidden
 * by 2. 4 and 6 are in the viewing plane while 5 juts out above plane and
 * 7 dips behind plane.
 *
 *  4   7 5 6
 *    : ·/ :
 *     (__)
 *    ·/  \·
 *   01    23
 *
 * Below is a top-down view. Strong lines indicate above-plane bonds, dots
 * indicate below-plane bonds.
 *
 *   0   7   3
 *     · | ·
 *   4 –( )– 6
 *     · | ·
 *   1   5   2
 *
 * @endverbatim
 *
 * Reference coodrinates:
 * - [W(CN)8]2-, (Alcock, M. W.; Samotus, A.; Szklarzewicz, J. J. Chem. Soc.,
 *   Dalton Trans. 1993, 885. DOI: 10.1039/dt9930000885, CSD: PECMOZ), took
 *   idealized geometry from http://symmetry.otterbein.edu/gallery/ instead
 *   of distorted crystallographic coordinates, normalized lengths
 *
 */
struct SquareAntiprism {
  static constexpr Shape shape = Shape::SquareAntiprism;
  static constexpr PointGroup pointGroup = PointGroup::D4d;
  static constexpr unsigned size = 8;
  static constexpr char stringName[] = "square antiprism";
  static constexpr std::array<temple::Vector, 8> coordinates {{
    // [W(CN)8]2-, idealized to square antiprism
    {-0.23838567,  0.50141283,  0.83171957},
    {-0.7568846,   0.61167543, -0.2301714 },
    { 0.3080136,   0.58106771, -0.75331795},
    { 0.82651172,  0.47080587,  0.30857773},
    {-0.79018301, -0.51909014,  0.32581627},
    {-0.39653401, -0.46341671, -0.79246813},
    { 0.72055552, -0.56338997, -0.40421711},
    { 0.32690564, -0.61906403,  0.71406753}
  }};

/*!
 * An upper triangular matrix containing angles between particules i,j in
 * degrees using the square antiprismatic reference coordinates
 */
  static constexpr auto angleLookupTable = temple::makeUpperTriangularMatrix(
    detail::makeArray<size>(coordinates)
  );

  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return angleLookupTable.at(
      std::min(a, b),
      std::max(a, b)
    );
  }
  static constexpr std::array<
    std::array<unsigned, 8>,
    2
  > rotations {{
    {{3, 0, 1, 2, 7, 4, 5, 6}}, // C4 axial
    /* 180° on equatorial axis in plane with 4, 6
     *
     *   1   5   2
     *     \ · /
     * – 4 ·( )· 6 –––––– equatorial axis
     *     / · \
     *   0   7   3
     *
     * and 45° anticlockwise on axial axis
     *
     *   5   2   6
     *     · | ·
     *   1 –( )– 3
     *     · | ·
     *   4   0   7
     *
     */
    {{5, 4, 7, 6, 1, 0, 3, 2}},
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
      {{0, 1, 4, 6}},
      {{1, 2, 5, 7}},
      {{2, 3, 6, 4}},
      {{3, 0, 7, 5}},
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    8
  > tetrahedra {{
    {{7, 0, 4, ORIGIN_PLACEHOLDER}},
    {{0, 4, ORIGIN_PLACEHOLDER, 1}},
    {{4, 1, 5, ORIGIN_PLACEHOLDER}},
    {{1, 5, ORIGIN_PLACEHOLDER, 2}},
    {{5, 2, 6, ORIGIN_PLACEHOLDER}},
    {{2, 6, ORIGIN_PLACEHOLDER, 3}},
    {{6, 3, 7, ORIGIN_PLACEHOLDER}},
    {{3, 7, ORIGIN_PLACEHOLDER, 0}}
  }};
#endif
  static constexpr std::array<unsigned, 8> mirror {{2, 1, 0, 3, 5, 4, 7, 6}};
};

/*!
 * @brief A cube shape
 */
struct Cube {
  static constexpr Symmetry::Shape shape = Symmetry::Shape::Cube;
  static constexpr unsigned size = 8;
  static constexpr char stringName[] = "cube";
  //! [V(CO)7]+ in C2v
  static constexpr std::array<temple::Vector, 8> coordinates {{
    {  0.577350,  0.577350,  0.577350},
    {  0.577350, -0.577350,  0.577350},
    {  0.577350, -0.577350, -0.577350},
    {  0.577350,  0.577350, -0.577350},
    { -0.577350,  0.577350,  0.577350},
    { -0.577350, -0.577350,  0.577350},
    { -0.577350, -0.577350, -0.577350},
    { -0.577350,  0.577350, -0.577350}
  }};
  static constexpr auto angleLookupTable = temple::makeUpperTriangularMatrix(
    detail::makeArray<size>(coordinates)
  );
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return angleLookupTable.at(
      std::min(a, b),
      std::max(a, b)
    );
  }
  static constexpr std::array<
    std::array<unsigned, 8>,
    2
  > rotations {{
    {{3, 0, 1, 2, 7, 4, 5, 6}}, // C4
    {{0, 4, 5, 1, 3, 7, 6, 2}} // C3
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, 1, 3, 5}},
    {{2, 4, 6, 7}}
  }};
  static constexpr std::array<unsigned, 8> mirror {{1, 0, 3, 2, 5, 4, 7, 6}};
};

/*!
 * @brief A bicapped trigonal prism
 */
struct BicappedTrigonalPrism {
  static constexpr Symmetry::Shape shape = Symmetry::Shape::BicappedTrigonalPrism;
  static constexpr unsigned size = 8;
  static constexpr char stringName[] = "bicapped trigonal prism";
  static constexpr std::array<temple::Vector, 8> coordinates {{
    {  0.577350,  0.577350,  0.577350},
    {  0.577350, -0.577350,  0.577350},
    {  0.577350, -0.577350, -0.577350},
    {  0.577350,  0.577350, -0.577350},
    { -0.577350,  0.577350,  0.577350},
    { -0.577350, -0.577350,  0.577350},
    { -0.577350, -0.577350, -0.577350},
    { -0.577350,  0.577350, -0.577350}
  }};
  static constexpr auto angleLookupTable = temple::makeUpperTriangularMatrix(
    detail::makeArray<size>(coordinates)
  );
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    return angleLookupTable.at(
      std::min(a, b),
      std::max(a, b)
    );
  }
  static constexpr std::array<
    std::array<unsigned, 8>,
    2
  > rotations {{
    {{3, 0, 1, 2, 7, 4, 5, 6}}, // C4
    {{0, 4, 5, 1, 3, 7, 6, 2}} // C3
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{0, 1, 3, 5}},
    {{2, 4, 6, 7}}
  }};
  static constexpr std::array<unsigned, 8> mirror {{1, 0, 3, 2, 5, 4, 7, 6}};
};

//! Type collecting all types of the Symmetry classes.
using allShapeDataTypes = std::tuple<
  Line,
  Bent,
  EquilateralTriangle, // 3
  ApicalTrigonalPyramid,
  T,
  Tetrahedron, // 4
  Square,
  Disphenoid,
  TrigonalPyramid,
  SquarePyramid, // 5
  TrigonalBipyramid,
  Pentagon,
  Octahedron, // 6
  TrigonalPrism,
  PentagonalPyramid,
  PentagonalBipyramid, // 7
  SquareAntiprism // 8
>;

static_assert(
  temple::TupleType::allOf<allShapeDataTypes, concepts::ShapeClass>(),
  "Not all shape data types fulfill the ShapeClass concept"
);

} // namespace data

} // namespace Symmetry

} // namespace Scine

#endif
