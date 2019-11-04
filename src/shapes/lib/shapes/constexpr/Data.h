/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Central symmetry data class definitions
 */

#ifndef INCLUDE_MOLASSEMBLER_SHAPES_CONSTEXPR_DATA_H
#define INCLUDE_MOLASSEMBLER_SHAPES_CONSTEXPR_DATA_H

#include "temple/constexpr/Vector.h"
#include "temple/constexpr/TupleType.h"

#include "shapes/Shapes.h"
#include "shapes/PointGroups.h"
#include "shapes/constexpr/CompileTimeOptions.h"
#include "shapes/constexpr/AngleLookup.h"

namespace Scine {

namespace Shapes {

//! A placeholder value for constexpr tetrahedra specification of origin
constexpr unsigned ORIGIN_PLACEHOLDER = std::numeric_limits<unsigned>::max();

/*!
 * @brief All symmetry data classes and some minor helper functions
 *
 * Each symmetry data class follows a concept seen in the .cpp file
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
 * @brief Bent symmetry at 107°
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
 * @brief Equilateral triangle shape (planar)
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
 * misrepresented, but all positions including the central vertex are in one
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
 * @brief Mono-vacant tetrahedron shape
 *
 * This symmetry is widely called trigonal pyramidal, but that name risks being
 * confused with a face-centered trigonal pyramid. This name should be preferred
 * for this geometry.
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
struct VacantTetrahedron {
  static constexpr Shape shape = Shape::VacantTetrahedron;
  static constexpr PointGroup pointGroup = PointGroup::C3v;
  static constexpr unsigned size = 3;
  static constexpr char stringName[] = "vacant tetrahedron";
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
 * @brief A regular tetrahedron shape
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

    return 2 * temple::Math::atan(M_SQRT2);
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
 *    (_) <- central vertex
 *    / \
 *   0   1
 *
 * @endverbatim
 *
 * Once again, angles are misrepresented, this is really a square.
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
 * @brief A seesaw shape
 *
 * More precisely, this has an angle of 120° between substituents 1 and 2. That
 * makes this shape both an equatorially mono-vacant trigonal bipyramid and an
 * edge-centered tetragonal disphenoid, if either of those help you imagine the
 * shape better.
 *
 * @verbatim
 *
 *   0 – (_) – 3
 *       / :
 *      1   2
 *
 * @endverbatim
 */
struct Seesaw {
  static constexpr Shape shape = Shape::Seesaw;
  static constexpr PointGroup pointGroup = PointGroup::C2v;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "seesaw";
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
 * @brief A face-centered trigonal pyramid shape = trig. pl. + an axial ligand
 *
 * A trigonal planar shape + an axial ligand, or alternatively an axially
 * mono-vacant trigonal bipyramid.
 *
 * @verbatim
 *
 * Viewed diagonally. The central vertex is ( ), 3 is apical, 0, 1, 2 and the
 * central vertex are in a plane.
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
 * @brief A square pyramid shape, the J1 solid (central position is square-face centered)
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
 * @brief A trigonal bipyramid shape, the J12 solid
 *
 * @verbatim
 *
 * Viewed from the top of the pyramid. The central vertex is ( ), 3 and 4
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
 * @brief A pentagon shape (planar)
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
 * @brief A regular octahedron
 *
 * @verbatim
 *
 * The central vertex is ( ), 4 and 5 are axial, the rest equatorial.
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
 * @brief A regular trigonal prism shape
 *
 * Squares and equilateral triangles as faces.
 *
 * @verbatim
 *
 *   2  0  1
 *    : | :
 *     (_)
 *    : | :
 *   5  3  3
 *
 * Where /, \ denote bonds in front of the view plane, : denotes bonds
 * behind the view plane.
 *
 * @endverbatim
 */
struct TrigonalPrism {
  static constexpr Shape shape = Shape::TrigonalPrism;
  static constexpr PointGroup pointGroup = PointGroup::D3h;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "trigonal prism";
  static constexpr std::array<temple::Vector, 6> coordinates {{
    { 0.755929,  0.000000,  0.654654},
    {-0.377964,  0.654654,  0.654654},
    {-0.377964, -0.654654,  0.654654},
    { 0.755929,  0.000000, -0.654654},
    {-0.377964,  0.654654, -0.654654},
    {-0.377964, -0.654654, -0.654654}
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
    std::array<unsigned, 6>,
    2
  > rotations {{
    {{2, 0, 1, 5, 3, 4}}, // C3 axial
    {{3, 5, 4, 0, 2, 1}} // C2 betw. 0, 3
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {{ORIGIN_PLACEHOLDER, 0, 2, 1}},
    {{3, ORIGIN_PLACEHOLDER, 5, 4}}
  }};
  static constexpr std::array<unsigned, 6> mirror {{0, 2, 1, 3, 5, 4}};
};

/*!
 * @brief A pentagonal pyramid shape, the J2 solid
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
 * ( ) signifies the central vertex beneath it
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
 * @brief Hexagon shape (planar)
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
 * @brief A pentagonal bipyramid shape, the J13 solid
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
 *
 * This is a gyroelongated triangular pyramid, or alternatively a "capped
 * triangular antiprism", depending on whatever helps you visualize it.
 */
struct CappedOctahedron {
  static constexpr Shape shape = Shape::CappedOctahedron;
  static constexpr PointGroup pointGroup = PointGroup::C3v;
  static constexpr unsigned size = 7;
  static constexpr char stringName[] = "capped octahedron";
  /*! Spherized [V(CO)7]+ in C3v (find local minimium of Thomson potential)
   *
   * from Jay W. Dicke, Nathan J. Stibrich, Henry F. Schaefer,
   * V(CO)7+: A capped octahedral structure completes the 18-electron rule,
   * Chemical Physics Letters, Volume 456, Issues 1–3, 2008.
   */
  static constexpr std::array<temple::Vector, 7> coordinates {{
    { 0.000000,  0.000000,  1.000000},
    { 0.957729,  0.000000,  0.287673},
    {-0.478864,  0.829418,  0.287673},
    {-0.478864, -0.829418,  0.287673},
    { 0.389831,  0.675207, -0.626200},
    {-0.779662,  0.000000, -0.626200},
    { 0.389831, -0.675207, -0.626200}
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
 * @brief A capped trigonal prism shape, spherized J49 solid in C2v
 *
 * Square-face capped.
 */
struct CappedTrigonalPrism {
  static constexpr Shape shape = Shape::CappedTrigonalPrism;
  static constexpr PointGroup pointGroup = PointGroup::C2v;
  static constexpr unsigned size = 7;
  static constexpr char stringName[] = "capped trigonal prism";
  /*! [V(CO)7]+ in C2v, from same source as CappedOctahedron
   *
   * Minimized to local minimum in Thomson potential
   */
  static constexpr std::array<temple::Vector, 7> coordinates {{
    { -0.000000, -0.000000,  1.000000},
    {  0.984798, -0.069552,  0.159173},
    { -0.069552,  0.984798,  0.159173},
    { -0.984798,  0.069552,  0.159173},
    {  0.069552, -0.984798,  0.159173},
    {  0.413726,  0.413726, -0.810964},
    { -0.413726, -0.413726, -0.810964}
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
 * @brief Regular square antiprism shape
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
 */
struct SquareAntiprism {
  static constexpr Shape shape = Shape::SquareAntiprism;
  static constexpr PointGroup pointGroup = PointGroup::D4d;
  static constexpr unsigned size = 8;
  static constexpr char stringName[] = "square antiprism";
  static constexpr std::array<temple::Vector, 8> coordinates {{
    { 0.607781,  0.607781,  0.511081},
    {-0.607781,  0.607781,  0.511081},
    {-0.607781, -0.607781,  0.511081},
    { 0.607781, -0.607781,  0.511081},
    { 0.859533,  0.000000, -0.511081},
    { 0.000000,  0.859533, -0.511081},
    {-0.859533,  0.000000, -0.511081},
    {-0.000000, -0.859533, -0.511081}
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

  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
      {{0, 1, 4, 6}},
      {{1, 2, 5, 7}},
      {{2, 3, 6, 4}},
      {{3, 0, 7, 5}},
  }};
  static constexpr std::array<unsigned, 8> mirror {{2, 1, 0, 3, 5, 4, 7, 6}};
};

/*!
 * @brief A regular cube
 */
struct Cube {
  static constexpr Shape shape = Shape::Cube;
  static constexpr PointGroup pointGroup = PointGroup::Oh;
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
    {{4, 5, 1, 0, 7, 6, 2, 3}} // C4'
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

/**
 * @brief Trigonal dodecahedron, snub disphenoid shape, spherized J84 solid in D2d
 */
struct TrigonalDodecahedron {
  static constexpr Shape shape = Shape::TrigonalDodecahedron;
  static constexpr PointGroup pointGroup = PointGroup::D2d;
  static constexpr unsigned size = 8;
  static constexpr char stringName[] = "trigonal dodecahedron";
  static constexpr std::array<temple::Vector, 8> coordinates {{
    {  0.620913,  0.000000, -0.783880},
    { -0.620913,  0.000000, -0.783880},
    {  0.000000,  0.620913,  0.783880},
    { -0.000000, -0.620913,  0.783880},
    {  0.950273,  0.000000,  0.311417},
    { -0.950273,  0.000000,  0.311417},
    {  0.000000,  0.950273, -0.311417},
    {  0.000000, -0.950273, -0.311417}
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
    {1, 0, 3, 2, 5, 4, 7, 6}, // C2z between 01
    {2, 3, 0, 1, 6, 7, 4, 5} // C2x + C4z
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {4, 2, 3, 5},
    {0, 6, 7, 1}
  }};
  static constexpr std::array<unsigned, 8> mirror {{0, 1, 3, 2, 4, 5, 7, 6}};
};

/**
 * @brief Hexagonal bipyramid shape
 *
 * Indices 0-5 ccw in the equatorial plane, 6 above, 7 below
 */
struct HexagonalBipyramid {
  static constexpr Shape shape = Shape::HexagonalBipyramid;
  static constexpr PointGroup pointGroup = PointGroup::D6h;
  static constexpr unsigned size = 8;
  static constexpr char stringName[] = "hexagonal bipyramid";
  static constexpr std::array<temple::Vector, 8> coordinates {{
    { 1.000000,  0.000000,  0.000000},
    { 0.500000,  0.866025,  0.000000},
    {-0.500000,  0.866025,  0.000000},
    {-1.000000,  0.000000,  0.000000},
    {-0.500000, -0.866025,  0.000000},
    { 0.500000, -0.866025,  0.000000},
    { 0.000000,  0.000000,  1.000000},
    { 0.000000,  0.000000, -1.000000}
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
    {5, 0, 1, 2, 3, 4, 6, 7}, // axial C6
    {0, 5, 4, 3, 2, 1, 7, 6} // C2 around 0-3
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > tetrahedra {{
    {6, 0, 1, 7},
    {6, 4, 5, 7},
    {6, 2, 3, 7}
  }};
  static constexpr std::array<unsigned, 8> mirror {{0, 5, 4, 3, 2, 1, 6, 7}};
};

/**
 * @brief Tricapped trigonal prism, spherized J51 solid in D3h
 *
 * Square-face tricapped. The coordinates are the solution to the Thomson
 * problem with 9 particles.
 */
struct TricappedTrigonalPrism {
  static constexpr Shape shape = Shape::TricappedTrigonalPrism;
  static constexpr PointGroup pointGroup = PointGroup::D3h;
  static constexpr unsigned size = 9;
  static constexpr char stringName[] = "tricapped trigonal prism";
  static constexpr std::array<temple::Vector, 9> coordinates {{
    { 0.914109572223, -0.182781178690, -0.361931942064},
    { 0.293329304506,  0.734642489361, -0.611766566546},
    {-0.480176899428, -0.046026929940,  0.875963279468},
    {-0.705684904851,  0.704780196051, -0.072757750931},
    { 0.370605109670,  0.769162968265,  0.520615194684},
    {-0.904030464226, -0.412626217894, -0.111662545460},
    {-0.162180419233, -0.247163999394, -0.955304908927},
    { 0.063327560246, -0.997971078243, -0.006583851785},
    { 0.610701141906, -0.322016246902,  0.723429092590}
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
    std::array<unsigned, 9>,
    2
  > rotations {{
    {7, 8, 3, 4, 2, 1, 0, 6, 5}, // C3 ccw between 2-4-3
    {2, 5, 0, 6, 7, 1, 3, 4, 8} // C2 at 8
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {1, 0, 6, 7},
    {5, 2, 3, 4}
  }};
  static constexpr std::array<unsigned, 9> mirror {{7, 5, 4, 3, 2, 1, 6, 0, 8}};
};

/**
 * @brief Capped square antiprism shape, spherized J10 solid in C4v
 *
 * Coordinates minimized to C4v local minimum of Thomson potential
 */
struct CappedSquareAntiprism {
  static constexpr Shape shape = Shape::CappedSquareAntiprism;
  static constexpr PointGroup pointGroup = PointGroup::C4v;
  static constexpr unsigned size = 9;
  static constexpr char stringName[] = "capped square antiprism";
  static constexpr std::array<temple::Vector, 9> coordinates {{
    { -0.000000,  0.932111,  0.362172},
    { -0.000000, -0.932111,  0.362172},
    {  0.932111, -0.000000,  0.362172},
    { -0.932111,  0.000000,  0.362172},
    {  0.559626,  0.559626, -0.611258},
    {  0.559626, -0.559626, -0.611258},
    { -0.559626,  0.559626, -0.611258},
    { -0.559626, -0.559626, -0.611258},
    {  0.000000,  0.000000,  1.000000}
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
    std::array<unsigned, 9>,
    1
  > rotations {{
    {2, 3, 1, 0, 5, 7, 4, 6, 8}
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    {6, 3, 0, 4},
    {7, 5, 1, 2}
  }};
  static constexpr std::array<unsigned, 9> mirror {{0, 1, 3, 2, 6, 7, 4, 5, 8}};
};

/**
 * @brief Capped square antiprism shape, spherized J10 solid in C4v
 */
struct HeptagonalBipyramid {
  static constexpr Shape shape = Shape::HeptagonalBipyramid;
  static constexpr PointGroup pointGroup = PointGroup::D7h;
  static constexpr unsigned size = 9;
  static constexpr char stringName[] = "heptagonal bipyramid";
  static constexpr std::array<temple::Vector, 9> coordinates {{
    { 1.000000,  0.000000,  0.000000},
    { 0.623490,  0.781831,  0.000000},
    {-0.222521,  0.974928,  0.000000},
    {-0.900969,  0.433884,  0.000000},
    {-0.900969, -0.433884,  0.000000},
    {-0.222521, -0.974928,  0.000000},
    { 0.623490, -0.781831,  0.000000},
    { 0.000000,  0.000000,  1.000000},
    { 0.000000,  0.000000, -1.000000}
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
    std::array<unsigned, 9>,
    2
  > rotations {{
    {6, 0, 1, 2, 3, 4, 5, 7, 8}, // axial C7
    {0, 6, 5, 4, 3, 2, 1, 8, 7} // C2 around 1 and between 4 and 5
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > tetrahedra {{
    {8, 1, 0, 7},
    {8, 3, 2, 7},
    {8, 5, 4, 7}
  }};
  static constexpr std::array<unsigned, 9> mirror {{0, 6, 5, 4, 3, 2, 1, 7, 8}};
};

/**
 * @brief Bicapped square antiprism shape, spherized J17 shape in D4h
 *
 * This is the solution to the Thomson problem with 10 particles.
 */
struct BicappedSquareAntiprism {
  static constexpr Shape shape = Shape::BicappedSquareAntiprism;
  static constexpr PointGroup pointGroup = PointGroup::D4h;
  static constexpr unsigned size = 10;
  static constexpr char stringName[] = "bicapped square antiprism";
  static constexpr std::array<temple::Vector, 10> coordinates {{
    { 0.978696890330,  0.074682616274,  0.191245663177},
    { 0.537258145625,  0.448413180814, -0.714338368164},
    {-0.227939324473, -0.303819959434, -0.925060590777},
    { 0.274577116268,  0.833436432027,  0.479573895237},
    {-0.599426405232,  0.240685139624,  0.763386303437},
    {-0.424664555168,  0.830194107787, -0.361161679833},
    {-0.402701180119, -0.893328907767,  0.199487398294},
    { 0.552788606831, -0.770301636525, -0.317899583084},
    { 0.290107593166, -0.385278374104,  0.876012647646},
    {-0.978696887344, -0.074682599351, -0.191245685067}
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
    std::array<unsigned, 10>,
    2
  > rotations {{
    {0, 7, 6, 1, 5, 2, 4, 8, 3, 9}, // C4z (0-9)
    {9, 5, 3, 2, 7, 1, 8, 4, 6, 0} // C2x + C8z
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > tetrahedra {{
    {0, 6, 8, 7},
    {9, 3, 4, 5},
    {9, 2, 1, 5}
  }};
  static constexpr std::array<unsigned, 10> mirror {{0, 1, 5, 7, 6, 2, 4, 3, 8, 9}};
};

/**
 * @brief Edge contracted icosahedron shape
 *
 * This is the solution for 11 particles in the Thomson problem in C2v point
 * group symmetry.
 */
struct EdgeContractedIcosahedron {
  static constexpr Shape shape = Shape::EdgeContractedIcosahedron;
  static constexpr PointGroup pointGroup = PointGroup::C2v;
  static constexpr unsigned size = 11;
  static constexpr char stringName[] = "edge-contracted icosahedron";
  static constexpr std::array<temple::Vector, 11> coordinates {{
    { 0.153486836562, -0.831354332797,  0.534127105044},
    { 0.092812115769,  0.691598091278, -0.716294626049},
    { 0.686120068086,  0.724987503180,  0.060269166267},
    { 0.101393837471,  0.257848797505,  0.960850293931},
    {-0.143059218646, -0.243142754178, -0.959382958495},
    {-0.909929380017,  0.200934944687, -0.362841110384},
    {-0.405338453688,  0.872713317547,  0.272162090194},
    { 0.896918545883, -0.184616420020,  0.401813264476},
    { 0.731466092268, -0.415052523977, -0.541007170195},
    {-0.439821168531, -0.864743799130, -0.242436592901},
    {-0.773718984882, -0.203685975092,  0.599892453681}
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
    std::array<unsigned, 11>,
    1
  > rotations {{
    {1, 0, 9, 5, 7, 3, 10, 4, 8, 2, 6} // C2
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > tetrahedra {{
    {6,  1, 5, 4},
    {3, 10, 0, 9},
    {1,  2, 8, 7}
  }};
  static constexpr std::array<unsigned, 11> mirror {{2, 9, 0, 3, 4, 5, 10, 7, 8, 1, 6}};
};

/**
 * @brief Regular icosahedron shape
 */
struct Icosahedron {
  static constexpr Shape shape = Shape::Icosahedron;
  static constexpr PointGroup pointGroup = PointGroup::Ih;
  static constexpr unsigned size = 12;
  static constexpr char stringName[] = "icosahedron";
  static constexpr std::array<temple::Vector, 12> coordinates {{
    { 0.525731,  0.000000,  0.850651},
    { 0.525731,  0.000000, -0.850651},
    {-0.525731,  0.000000,  0.850651},
    {-0.525731,  0.000000, -0.850651},
    { 0.850651,  0.525731,  0.000000},
    { 0.850651, -0.525731,  0.000000},
    {-0.850651,  0.525731,  0.000000},
    {-0.850651, -0.525731,  0.000000},
    { 0.000000,  0.850651,  0.525731},
    { 0.000000,  0.850651, -0.525731},
    { 0.000000, -0.850651,  0.525731},
    { 0.000000, -0.850651, -0.525731}
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
    std::array<unsigned, 12>,
    3
  > rotations {{
    {0, 11, 8, 3, 5, 10, 9, 6, 4, 1, 2, 7}, // C5 around 0-3
    {8, 5, 6, 11, 4, 0, 3, 7, 9, 1, 2, 10}, // C5 around 4-7
    {2, 3, 0, 1, 7, 6, 5, 4, 10, 11, 8, 9} // C2 between 0-2 / 1-3
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
    {0, 2, 10, 8},
    {1, 3, 9, 11},
    {1, 4, 5, 0},
    {3, 7, 6, 2}
  }};
  static constexpr std::array<unsigned, 12> mirror {{0, 1, 2, 3, 5, 4, 7, 6, 10, 11, 8, 9}};
};

/**
 * @brief Regular cuboctahedron shape with Oh symmetry
 */
struct Cuboctahedron {
  static constexpr Shape shape = Shape::Cuboctahedron;
  static constexpr PointGroup pointGroup = PointGroup::Oh;
  static constexpr unsigned size = 12;
  static constexpr char stringName[] = "cuboctahedron";
  static constexpr std::array<temple::Vector, 12> coordinates {{
    { 0.707107,  0.000000,  0.707107},
    { 0.707107,  0.000000, -0.707107},
    {-0.707107,  0.000000,  0.707107},
    {-0.707107,  0.000000, -0.707107},
    { 0.707107,  0.707107,  0.000000},
    { 0.707107, -0.707107,  0.000000},
    {-0.707107,  0.707107,  0.000000},
    {-0.707107, -0.707107,  0.000000},
    { 0.000000,  0.707107,  0.707107},
    { 0.000000,  0.707107, -0.707107},
    { 0.000000, -0.707107,  0.707107},
    { 0.000000, -0.707107, -0.707107}
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
    std::array<unsigned, 12>,
    3
  > rotations {{
    {10, 11, 8, 9, 5, 7, 4, 6, 0, 1, 2, 3}, // C4 ccw 0-8-2-10
    {2, 0, 3, 1, 8, 10, 9, 11, 6, 4, 7, 5}, // C4 ccw 4-9-6-8
    {7, 6, 5, 4, 3, 2, 1, 0, 11, 9, 10, 8}, // C2 along 9-10
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > tetrahedra {{
    {ORIGIN_PLACEHOLDER, 6,  9,  8},
    {ORIGIN_PLACEHOLDER, 4,  1,  0},
    {ORIGIN_PLACEHOLDER, 5, 11, 10},
    {ORIGIN_PLACEHOLDER, 7,  3,  2}
  }};
  static constexpr std::array<unsigned, 12> mirror {{8, 9, 10, 11, 4, 6, 5, 7, 0, 1, 2, 3}};
};

//! Type collecting all types of the Symmetry classes.
using allShapeDataTypes = std::tuple<
  Line,
  Bent,
  EquilateralTriangle, // 3
  VacantTetrahedron,
  T,
  Tetrahedron, // 4
  Square,
  Seesaw,
  TrigonalPyramid,
  SquarePyramid, // 5
  TrigonalBipyramid,
  Pentagon,
  Octahedron, // 6
  TrigonalPrism,
  PentagonalPyramid,
  Hexagon,
  PentagonalBipyramid, // 7
  CappedOctahedron,
  CappedTrigonalPrism,
  SquareAntiprism, // 8
  Cube,
  TrigonalDodecahedron,
  HexagonalBipyramid,
  TricappedTrigonalPrism, // 9
  CappedSquareAntiprism,
  HeptagonalBipyramid,
  BicappedSquareAntiprism, // 10
  EdgeContractedIcosahedron, // 11
  Icosahedron, // 12,
  Cuboctahedron
>;

} // namespace data

} // namespace Shapes

} // namespace Scine

#endif
