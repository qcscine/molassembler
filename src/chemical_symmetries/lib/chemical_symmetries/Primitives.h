/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Central symmetry data class definitions
 */

#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_PRIMITIVES_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_PRIMITIVES_H

#include "temple/constexpr/Vector.h"

#include "chemical_symmetries/Names.h"
#include "chemical_symmetries/CompileTimeOptions.h"

#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
#include "chemical_symmetries/AngleLookup.h"
#endif

#include <cmath>

namespace Scine {

namespace Symmetry {

//! A placeholder value for constexpr tetrahedra specification of origin
constexpr unsigned ORIGIN_PLACEHOLDER = std::numeric_limits<unsigned>::max();

/*!
 * Namespace containing all symmetry data classes and some minor helper functions
 *
 * Each symmetry data class must have the following members, all of which must
 * be static constexpr (or static const in the exception of stringName):
 * - name (Symmetry::Name)
 * - size (unsigned)
 * - stringName (string)
 * - angleFunction ( double(const unsigned, const unsigned) )
 * - coordinates ( array<temple::Vector, N> ): N is size of symmetry,
 *   see above
 * - rotations ( array< array<unsigned, size>, R> ): R is however many
 *   rotations are needed to produce all superimposable rotations
 * - tetrahedra ( array< array<unsigned, 4>, T> ): T is however many tetrahedra
 *   are required to completely define the chirality
 * - mirror (array<unsigned, N>): plane may not slice molecule in any way, it
 *   is auxiliary to finding enantiomeric pairs. If a mirror were to yield a
 *   superposable stereopermutation, represent it by an empty array.
 *
 */
namespace data {

/*!
 * @brief Linear symmetry
 *
 * @verbatim
 *
 *  0 – (_) – 1
 *
 * @endverbatim
 */
struct Linear {
  static constexpr Symmetry::Name name = Symmetry::Name::Linear;
  static constexpr unsigned size = 2;
  static constexpr char stringName[] = "linear";
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
  static constexpr Symmetry::Name name = Symmetry::Name::Bent;
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
struct TrigonalPlanar {
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPlanar;
  static constexpr unsigned size = 3;
  static constexpr char stringName[] = "trigonal planar";
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
 * @brief A Tetrahedral symmetry missing a ligand
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
struct CutTetrahedral {
  static constexpr Symmetry::Name name = Symmetry::Name::CutTetrahedral;
  static constexpr unsigned size = 3;
  static constexpr char stringName[] = "cut tetrahedral";
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
struct TShaped {
  static constexpr Symmetry::Name name = Symmetry::Name::TShaped;
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
struct Tetrahedral {
  static constexpr Symmetry::Name name = Symmetry::Name::Tetrahedral;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "tetrahedral";
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
struct SquarePlanar {
  static constexpr Symmetry::Name name = Symmetry::Name::SquarePlanar;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "square planar";
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
 * @brief A seesaw symmetry
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
  static constexpr Symmetry::Name name = Symmetry::Name::Seesaw;
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
struct TrigonalPyramidal {
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPyramidal;
  static constexpr unsigned size = 4;
  static constexpr char stringName[] = "trigonal pyramidal";
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
struct SquarePyramidal {
  static constexpr Symmetry::Name name = Symmetry::Name::SquarePyramidal;
  static constexpr unsigned size = 5;
  static constexpr char stringName[] = "square pyramidal";
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
struct TrigonalBiPyramidal {
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalBiPyramidal;
  static constexpr unsigned size = 5;
  static constexpr char stringName[] = "trigonal bipyramidal";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    unsigned smaller = std::min(a, b), larger = std::max(a, b);
    if(larger < 3) {
      // -> smaller < 2, this means either 0,1 0,2 1,2 axial
      return temple::Math::toRadians<double>(120);
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
 * @brief A pentagonal planar symmetry
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
struct PentagonalPlanar {
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalPlanar;
  static constexpr unsigned size = 5;
  static constexpr char stringName[] = "pentagonal planar";
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
struct Octahedral {
  static constexpr Symmetry::Name name = Symmetry::Name::Octahedral;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "octahedral";
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
struct TrigonalPrismatic {
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPrismatic;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "trigonal prismatic";
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
struct PentagonalPyramidal {
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalPyramidal;
  static constexpr unsigned size = 6;
  static constexpr char stringName[] = "pentagonal pyramidal";
  static constexpr double angleFunction(const unsigned a, const unsigned b) {
    if(a == b) {
      return 0;
    }

    if(a == 5 || b == 5) {
      return M_PI / 2;
    }

    // remainder are identical to PentagonalPlanar
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
struct PentagonalBiPyramidal {
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalBiPyramidal;
  static constexpr unsigned size = 7;
  static constexpr char stringName[] = "pentagonal bipyramidal";
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

    // remainder are equatorial angles, like PentagonalPlanar
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
struct SquareAntiPrismatic {
  static constexpr Symmetry::Name name = Symmetry::Name::SquareAntiPrismatic;
  static constexpr unsigned size = 8;
  static constexpr char stringName[] = "square antiprismatic";
  static constexpr std::array<temple::Vector, 8> coordinates {{
    // [W(CN)8]2-, idealized to square antiprism
    {-0.23838567,  0.50141283,  0.83171957},
    {-0.7568846,   0.61167543, -0.2301714 },
    { 0.3080136,   0.58106771, -0.75331795},
    { 0.82651172,  0.47080587,  0.30857773},
    {-0.79018301, -0.51909014,  0.32581627},
    {-0.39653401, -0.46341671, -0.79246813},
    { 0.72055552, -0.56338997, -0.40421711},
    { 0.32690564, -0.61906403,  0.71406753},
  }};

#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
/*!
 * An upper triangular matrix containing angles between particules i,j in
 * degrees using the square antiprismatic reference coordinates
 */
  static constexpr auto angleLookupTable = temple::makeUpperTriangularMatrix(
    detail::makeArray<size>(coordinates)
  );
#endif

  static constexpr double angleFunction(const unsigned a, const unsigned b) {
#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
    if(a == b) {
      return 0;
    }

    return angleLookupTable.at(
      std::min(a, b),
      std::max(a, b)
    );
#else
    if(a == b) {
      return 0;
    }

    if(
      (a < 4 && b < 4)
      || (a >= 4 && b >= 4)
    ) { // in plane
      if((a + b) % 2 == 1) { // cis
        return temple::Math::toRadians<double>(72.9875);
      }

      // otherwise trans
      return temple::Math::toRadians<double>(114.475);
    }

    // remaining cases are between planes
    unsigned minDiff = std::min(a - b, b - a);
    if(minDiff == 3 || minDiff == 4 || minDiff == 7) { // short
      return temple::Math::toRadians<double>(78.05);
    }

    // last case is long between planes
    return temple::Math::toRadians<double>(142.275);
#endif
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


//! Type collecting all types of the Symmetry classes.
using allSymmetryDataTypes = std::tuple<
  Linear, // 2
  Bent,
  TrigonalPlanar, // 3
  CutTetrahedral,
  TShaped,
  Tetrahedral, // 4
  SquarePlanar,
  Seesaw,
  TrigonalPyramidal,
  SquarePyramidal, // 5
  TrigonalBiPyramidal,
  PentagonalPlanar,
  Octahedral, // 6
  TrigonalPrismatic,
  PentagonalPyramidal,
  PentagonalBiPyramidal, // 7
  SquareAntiPrismatic // 8
>;

} // namespace data

} // namespace Symmetry

} // namespace Scine

#endif
