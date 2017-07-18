#ifndef LIB_INCLUDE_SYMMETRIES_H
#define LIB_INCLUDE_SYMMETRIES_H

/* If USE_ALTERNATE_TETRAHEDRA is defined, a reduced set of tetrahedra
 * is used to subdivide higher symmetries. This may provide less information 
 * about the geometry when used but should improve performance as fewer 
 * tetrahedron volumes must be calculated.
 */
//#define USE_ALTERNATE_TETRAHEDRA

/* If USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE is defined, a table of all
 * angles resulting from a predefined set of positions is generated and that
 * symmetry's angle function turns into what is essentially a lookup table.
 */
#define USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE

#include "Eigen/Core"
#include "boost/optional.hpp"

#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
#include "ConstexprAngles.h"
#endif

#include "constexpr_magic/Math.h"
#include "template_magic/Containers.h"

#include <map>
#include <vector>
#include <functional>
#include <algorithm>

/* TODO
 * - Debug and Release builds
 * - Improve trigonal pyramidal coordinates definition to get 107.5 angle as a
 *   parameter.  Currently, the rotation angle choice of 111.5 works well, but
 *   completely arbitrary!
 * - Consider making constexpr calculation of all angles from coordinates into
 *   const lookup table
 * - Could replicate angle parametrization of coordinates with a constexpr
 *   matrix class and matrix * vector multiplication
 */

namespace Symmetry {

/* Typedefs */
using RotationsList = std::vector<
  std::vector<unsigned>
>;

/* All angle functions can be called with arbitrary (valid) parameters
 * without failing. Valid here means that a != b and less than the size of
 * the symmetry requested.
 *
 * They return angles in radians.
 */
using AngleFunctionType = std::function<
  double(const unsigned&, const unsigned&)
>;

/* All symmetries have a guess implementation of what could work as the defined
 * tetrahedra. Have to use boost::none to signal to replace this position with 
 * the central atom as it is not part of the indexing scheme used here.
 *
 * In case all higher symmetries than trigonal pyramidal are representable 
 * without boost::none and that proves to work, then perhaps make an exception 
 * for it and treat all others without the optional. If that cannot be done, 
 * consider refactoring (changing the numbering scheme in some fashion that 
 * boost::none does not have to be used.
 */
using TetrahedronList = std::vector<
  std::array<
    boost::optional<unsigned>,
    4
  >
>;

using CoordinateList = std::vector<Eigen::Vector3d>;

// Dynamic symmetry information class
struct SymmetryInformation {
  const std::string stringName;
  const unsigned size;
  const RotationsList rotations;
  const TetrahedronList tetrahedra;
  const CoordinateList coordinates;

  // Direct initialization
  SymmetryInformation(
    std::string stringName,
    unsigned size,
    RotationsList rotations,
    TetrahedronList tetrahedra,
    CoordinateList coordinates
  ) : stringName(stringName),
      size(size),
      rotations(rotations),
      tetrahedra(tetrahedra),
      coordinates(coordinates)
  {}
};

// Symmetry names list
enum class Name : unsigned {
  Linear, // 2
  Bent,
  TrigonalPlanar, // 3
  TrigonalPyramidal,
  TShaped,
  Tetrahedral, // 4
  SquarePlanar,
  Seesaw,
  SquarePyramidal, // 5
  TrigonalBiPyramidal,
  PentagonalPlanar,
  Octahedral, // 6
  TrigonalPrismatic,
  PentagonalPyramidal,
  PentagonalBiPyramidal, // 7
  SquareAntiPrismatic // 8
};

// DATA
constexpr unsigned nSymmetries = 16;

constexpr std::array<Name, nSymmetries> allNames {
  Name::Linear, // 2
  Name::Bent,
  Name::TrigonalPlanar, // 3
  Name::TrigonalPyramidal,
  Name::TShaped,
  Name::Tetrahedral, // 4
  Name::SquarePlanar,
  Name::Seesaw,
  Name::SquarePyramidal, // 5
  Name::TrigonalBiPyramidal,
  Name::PentagonalPlanar,
  Name::Octahedral, // 6
  Name::TrigonalPrismatic,
  Name::PentagonalPyramidal,
  Name::PentagonalBiPyramidal, // 7
  Name::SquareAntiPrismatic // 8
};

constexpr unsigned replaceMe = std::numeric_limits<unsigned>::max();

namespace data {

struct Linear {
  static constexpr Symmetry::Name name = Symmetry::Name::Linear;
  static constexpr unsigned size = 2;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return M_PI;
  }
  static constexpr std::array<ConstexprMagic::Vector, 2> coordinates {{
    { 1 , 0, 0 },
    { -1, 0, 0 }
  }};
  static constexpr std::array<
    std::array<unsigned, 2>,
    1
  > rotations {{
    {1, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

constexpr decltype(Linear::name) Linear::name;
constexpr decltype(Linear::size) Linear::size;
constexpr decltype(Linear::coordinates) Linear::coordinates;
constexpr decltype(Linear::rotations) Linear::rotations;
constexpr decltype(Linear::tetrahedra) Linear::tetrahedra;


struct Bent {
  /*
   *  1
   *   \
   *    (_) – 0
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Bent;
  static constexpr unsigned size = 2;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    /* subject to a lot of variation, between 90 and 109 degrees pursuant to 
     * english wikipedia, using experimental data here to improve instances
     * of this geometry on e.g. O center would a big improvement to DG runs
     */
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(107); 
  }
  static constexpr std::array<ConstexprMagic::Vector, 2> coordinates {{
    {1., 0., 0.},
    {-0.292372, 0.956305, 0.}
  }};
  static constexpr std::array<
    std::array<unsigned, 2>,
    1
  > rotations {{
    {1, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

constexpr decltype(Bent::name) Bent::name;
constexpr decltype(Bent::size) Bent::size;
constexpr decltype(Bent::coordinates) Bent::coordinates;
constexpr decltype(Bent::rotations) Bent::rotations;
constexpr decltype(Bent::tetrahedra) Bent::tetrahedra;


struct TrigonalPlanar {
  /*
   *     0
   *     |
   *    (_)
   *   /   \
   *  1     2
   *
   * This is not quite ideal since the angles are thoroughly misrepresented, 
   * but all positions including the central atom are in one plane. The angles
   * are idealized as 120°.
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPlanar;
  static constexpr unsigned size = 3;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(120);
  }
  static constexpr std::array<ConstexprMagic::Vector, 3> coordinates {{
      {1, 0, 0},
      {-0.5, 0.866025, 0},
      {-0.5, -0.866025, 0}
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    2
  > rotations {{
    {1, 2, 0}, // C3
    {0, 2, 1} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

constexpr decltype(TrigonalPlanar::name) TrigonalPlanar::name;
constexpr decltype(TrigonalPlanar::size) TrigonalPlanar::size;
constexpr decltype(TrigonalPlanar::coordinates) TrigonalPlanar::coordinates;
constexpr decltype(TrigonalPlanar::rotations) TrigonalPlanar::rotations;
constexpr decltype(TrigonalPlanar::tetrahedra) TrigonalPlanar::tetrahedra;


struct TrigonalPyramidal {
  /*
   *     0
   *     |
   *    (_)
   *   /   \
   *  1     2
   *
   * This is not quite ideal since the angles are thoroughly misrepresented, 
   * but all positions including the central atom are in one plane. The angles
   * are idealized as 120°.
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPyramidal;
  static constexpr unsigned size = 3;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(107.5);
  }
  static constexpr std::array<ConstexprMagic::Vector, 3> coordinates {{
    {0, -0.366501, 0.930418},
    {0.805765, -0.366501, -0.465209},
    {-0.805765, -0.366501, -0.465209}
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    1
  > rotations {{
    {2, 0, 1} // C3
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{replaceMe, 0, 1, 2}}
  }};
};

constexpr decltype(TrigonalPyramidal::name) TrigonalPyramidal::name;
constexpr decltype(TrigonalPyramidal::size) TrigonalPyramidal::size;
constexpr decltype(TrigonalPyramidal::coordinates) TrigonalPyramidal::coordinates;
constexpr decltype(TrigonalPyramidal::rotations) TrigonalPyramidal::rotations;
constexpr decltype(TrigonalPyramidal::tetrahedra) TrigonalPyramidal::tetrahedra;


struct TShaped {
  /*
   * 0 – (_) – 2
   *      |
   *      1
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TShaped;
  static constexpr unsigned size = 3;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    if((a + b) % 2 == 1) {
      return M_PI / 2;
    } 

    return M_PI;
  }
  static constexpr std::array<ConstexprMagic::Vector, 3> coordinates {{
    {-1, -0, -0},
    {0, 1, 0},
    {1, 0, 0},
  }};
  static constexpr std::array<
    std::array<unsigned, 3>,
    1
  > rotations {{
    {2, 1, 0} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

constexpr decltype(TShaped::name) TShaped::name;
constexpr decltype(TShaped::size) TShaped::size;
constexpr decltype(TShaped::coordinates) TShaped::coordinates;
constexpr decltype(TShaped::rotations) TShaped::rotations;
constexpr decltype(TShaped::tetrahedra) TShaped::tetrahedra;


struct Tetrahedral {
  /* 
   *    1
   *    |
   *   (0) (0 is on top, ( ) signifies the central atom  beneath it
   *  /   \
   * 2     3
   *
   * Remember Newman projections? This is sort of supposed to be that.
   *
   * Alternatively:
   *
   *    0
   *    |
   *   (_)   
   *  /  \ °3
   * 1    2 
   *
   * Where /, \ denote in front of plane bonds, ° a behind the plane bond.
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Tetrahedral;
  static constexpr unsigned size = 4;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    return ConstexprMagic::Math::toRadians<double>(109.5);
  }
  static constexpr std::array<ConstexprMagic::Vector, 4> coordinates {{
    {0, 1, 0},
    {0, -0.333807, 0.942641},
    {0.816351, -0.333807, -0.471321},
    {-0.816351, -0.333807, -0.471321}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    4
  > rotations {{
    {0, 3, 1, 2}, // C4, 1
    {2, 1, 3, 0}, // C4, 2
    {3, 0, 2, 1}, // C4, 3
    {1, 2, 0, 3}  // C4, 4
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > tetrahedra {{
    {{0, 1, 2, 3}}
  }};
};

constexpr decltype(Tetrahedral::name) Tetrahedral::name;
constexpr decltype(Tetrahedral::size) Tetrahedral::size;
constexpr decltype(Tetrahedral::coordinates) Tetrahedral::coordinates;
constexpr decltype(Tetrahedral::rotations) Tetrahedral::rotations;
constexpr decltype(Tetrahedral::tetrahedra) Tetrahedral::tetrahedra;


struct SquarePlanar {
  /* 
   * 3   2
   *  \_/
   *  (_) <- central atom
   *  / \
   * 0   1
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::SquarePlanar;
  static constexpr unsigned size = 4;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
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
  static constexpr std::array<ConstexprMagic::Vector, 4> coordinates {{
    {1, 0, 0},
    {0, 1, 0},
    {-1, -0, -0},
    {-0, -1, -0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > rotations {{
    {3, 0, 1, 2}, // C4
    {1, 0, 3, 2}, // C2
    {3, 2, 1, 0}  // C2'
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

constexpr decltype(SquarePlanar::name) SquarePlanar::name;
constexpr decltype(SquarePlanar::size) SquarePlanar::size;
constexpr decltype(SquarePlanar::coordinates) SquarePlanar::coordinates;
constexpr decltype(SquarePlanar::rotations) SquarePlanar::rotations;
constexpr decltype(SquarePlanar::tetrahedra) SquarePlanar::tetrahedra;


struct Seesaw {
  /*
   * 0 – (_) – 3
   *     / :
   *    1   2
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Seesaw;
  static constexpr unsigned size = 4;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    const auto& smaller = std::min(a, b);
    const auto& larger = std::max(a, b);
    if(smaller == 0 && larger == 3) {
      return M_PI;
    }

    if(smaller == 1 && larger == 2) {
      return ConstexprMagic::Math::toRadians<double>(120);
    }

    return M_PI / 2;
  }
  static constexpr std::array<ConstexprMagic::Vector, 4> coordinates {{
    {0, 1, 0},
    {1, 0, 0},
    {-0.5, 0, -0.866025},
    {-0, -1, -0}
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    1
  > rotations {{
    {3, 2, 1, 0} // C2
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
    {{0, replaceMe, 1, 2}},
    {{replaceMe, 3, 1, 2}},
  }};
#endif
};

constexpr decltype(Seesaw::name) Seesaw::name;
constexpr decltype(Seesaw::size) Seesaw::size;
constexpr decltype(Seesaw::coordinates) Seesaw::coordinates;
constexpr decltype(Seesaw::rotations) Seesaw::rotations;
constexpr decltype(Seesaw::tetrahedra) Seesaw::tetrahedra;


struct SquarePyramidal {
  /* 
   * 3   2
   *  \_/
   *  (4)   
   *  / \
   * 0   1
   *
   * Viewed from the top of the pyramid. The central atom is ( ), 5 is axial.
   *
   * Alternatively,
   *
   *    4
   * 3  |  2
   *  : | :    <- behind view plane
   *   (_)
   *  /   \    <- in front of view plane
   * 0     1
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::SquarePyramidal;
  static constexpr unsigned size = 5;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
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
  static constexpr std::array<ConstexprMagic::Vector, 5> coordinates {{
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
    {3, 0, 1, 2, 4} // C4
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
    {{0, 1, 4, replaceMe}},
    {{1, 2, 4, replaceMe}},
    {{2, 3, 4, replaceMe}},
    {{3, 0, 4, replaceMe}}
  }};
#endif
};

constexpr decltype(SquarePyramidal::name) SquarePyramidal::name;
constexpr decltype(SquarePyramidal::size) SquarePyramidal::size;
constexpr decltype(SquarePyramidal::coordinates) SquarePyramidal::coordinates;
constexpr decltype(SquarePyramidal::rotations) SquarePyramidal::rotations;
constexpr decltype(SquarePyramidal::tetrahedra) SquarePyramidal::tetrahedra;


struct TrigonalBiPyramidal {
  /* Viewed from the top of the pyramid. The central atom is ( ), 3 and 4 
   * are axial.
   *
   *     3
   *     |  2
   *     | :    <- behind view plane
   * 0--(_)
   *     | \    <- in front of view plane
   *     |  1
   *     4
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalBiPyramidal;
  static constexpr unsigned size = 5;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    unsigned smaller = std::min(a, b), larger = std::max(a, b);
    if(larger < 3) {
      // -> smaller < 2, this means either 0,1 0,2 1,2 axial
      return ConstexprMagic::Math::toRadians<double>(120);
    } else if(larger == 3) {
      // -> smaller < 3, this means {1,2,3}, 3 
      return M_PI / 2;
    } else if(smaller < 3) {
      // now, larger must be 4 (process of elimination), so if a is not 3:
      return M_PI / 2;
    } else {
      // only case left: 3,4
      return M_PI;
    }
  }
  static constexpr std::array<ConstexprMagic::Vector, 5> coordinates {{
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
    {2, 0, 1, 3, 4}, // C3
    {0, 2, 1, 4, 3}, // C2 on 0
    {2, 1, 0, 4, 3}, // C2 on 1
    {1, 0, 2, 4, 3} // C2 on 2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
      {{0, 1, 3, 2}},
      {{0, 1, 2, 4}}
  }};
};

constexpr decltype(TrigonalBiPyramidal::name) TrigonalBiPyramidal::name;
constexpr decltype(TrigonalBiPyramidal::size) TrigonalBiPyramidal::size;
constexpr decltype(TrigonalBiPyramidal::coordinates) TrigonalBiPyramidal::coordinates;
constexpr decltype(TrigonalBiPyramidal::rotations) TrigonalBiPyramidal::rotations;
constexpr decltype(TrigonalBiPyramidal::tetrahedra) TrigonalBiPyramidal::tetrahedra;


struct PentagonalPlanar {
  /* 
   * All in plane:
   *
   *      0
   *  1.  |  .4
   *    °(_)°
   *    /   \
   *   2     3
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalPlanar;
  static constexpr unsigned size = 5;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * ConstexprMagic::Math::toRadians<double>(72);
  }
  static constexpr std::array<ConstexprMagic::Vector, 5> coordinates {{
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
    {4, 0, 1, 2, 3}, // C5
    {0, 4, 3, 2, 1} // C2
  }};
  static constexpr std::array<
    std::array<unsigned, 4>,
    0
  > tetrahedra {{}};
};

constexpr decltype(PentagonalPlanar::name) PentagonalPlanar::name;
constexpr decltype(PentagonalPlanar::size) PentagonalPlanar::size;
constexpr decltype(PentagonalPlanar::coordinates) PentagonalPlanar::coordinates;
constexpr decltype(PentagonalPlanar::rotations) PentagonalPlanar::rotations;
constexpr decltype(PentagonalPlanar::tetrahedra) PentagonalPlanar::tetrahedra;


struct Octahedral {
  /* The central atom is ( ), 4 and 5 are axial, the rest equatorial.
   *
   *     4
   *  3  |  2
   *   : | :
   *    (_)        
   *   / | \
   *  0  |  1
   *     5
   *
   * Where /, \ denote bonds in front of the view plane, : denotes bonds
   * behind the view plane.
   */
  static constexpr Symmetry::Name name = Symmetry::Name::Octahedral;
  static constexpr unsigned size = 6;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
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
  static constexpr std::array<ConstexprMagic::Vector, 6> coordinates {{
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
    {3, 0, 1, 2, 4, 5}, // vertical C4
    {0, 5, 2, 4, 1, 3}, // horizontal C4
    {4, 1, 5, 3, 2, 0} // horizontal C4'
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
    {{3, 0, 4, replaceMe}},
    {{0, 1, 4, replaceMe}},
    {{1, 2, 4, replaceMe}},
    {{2, 3, 4, replaceMe}},
    {{3, 0, replaceMe, 5}},
    {{0, 1, replaceMe, 5}},
    {{1, 2, replaceMe, 5}},
    {{2, 3, replaceMe, 5}}
  }};
#endif
};

constexpr decltype(Octahedral::name) Octahedral::name;
constexpr decltype(Octahedral::size) Octahedral::size;
constexpr decltype(Octahedral::coordinates) Octahedral::coordinates;
constexpr decltype(Octahedral::rotations) Octahedral::rotations;
constexpr decltype(Octahedral::tetrahedra) Octahedral::tetrahedra;


struct TrigonalPrismatic {
  /* 
   *  3  4  5
   *   : | :
   *    (_)        
   *   : | :
   *  0  1  2
   *
   * Where /, \ denote bonds in front of the view plane, : denotes bonds
   * behind the view plane.
   *
   * Angles 
   *  0-1, 0-2 -> 86°
   *  0-3 -> 76°
   *  0-4, 0-5 -> 134°
   * 
   * From [W(CH3)6], Haaland, Hammel, Rypdal, Volden, J. Am. Chem. Soc. 1990 
   */
  static constexpr Symmetry::Name name = Symmetry::Name::TrigonalPrismatic;
  static constexpr unsigned size = 6;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    } 
    
    // Between plane symmetric
    if(std::min(a - b, b - a) == 3) {
      return ConstexprMagic::Math::toRadians<double>(76);
    } 

    // In plane triangle
    if(
      (a < 3 && b < 3)
      || (a >= 3 && b >= 3)
    ) {
      return ConstexprMagic::Math::toRadians<double>(86);
    } 

    // Between plane asymmetric
    return ConstexprMagic::Math::toRadians<double>(134);
  }
  static constexpr std::array<ConstexprMagic::Vector, 6> coordinates {{
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
    {2, 0, 1, 5, 3, 4}, // C3 axial
    {5, 4, 3, 2, 1, 0} // C2 betw. 1, 4
  }};

  static constexpr std::array<
    std::array<unsigned, 4>,
    2
  > tetrahedra {{
    // TODO dubious if this captures all relevant information, too limited
    {{replaceMe, 0, 1, 2}},
    {{3, replaceMe, 4, 5}}
  }};
};

constexpr decltype(TrigonalPrismatic::name) TrigonalPrismatic::name;
constexpr decltype(TrigonalPrismatic::size) TrigonalPrismatic::size;
constexpr decltype(TrigonalPrismatic::coordinates) TrigonalPrismatic::coordinates;
constexpr decltype(TrigonalPrismatic::rotations) TrigonalPrismatic::rotations;
constexpr decltype(TrigonalPrismatic::tetrahedra) TrigonalPrismatic::tetrahedra;


struct PentagonalPyramidal {
  /* 
   *      0
   *  1.  |  .4
   *    °(5)°
   *    :   :
   *   2     3
   *
   * 0-4 are in plane,
   * 5 is above plane,
   * ( ) signifies the central atom beneath it
   */
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalPyramidal;
  static constexpr unsigned size = 6;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
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
    ) * ConstexprMagic::Math::toRadians<double>(72);
  }
  static constexpr std::array<ConstexprMagic::Vector, 6> coordinates {{
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
    {4, 0, 1, 2, 3, 5} // C5 axial
  }};

#ifdef USE_ALTERNATE_TETRAHEDRA
  static constexpr std::array<
    std::array<unsigned, 4>,
    3
  > tetrahedra {{
    {{0, 1, 5, 2}},
    {{2, 3, 5, 4}},
    {{4, 5, replaceMe, 0}}
  }};
#else // Regular
  static constexpr std::array<
    std::array<unsigned, 4>,
    5
  > tetrahedra {{
    {{0, replaceMe, 1, 5}},
    {{1, replaceMe, 2, 5}},
    {{2, replaceMe, 3, 5}},
    {{3, replaceMe, 4, 5}},
    {{4, replaceMe, 0, 5}}
  }};
#endif
};

constexpr decltype(PentagonalPyramidal::name) PentagonalPyramidal::name;
constexpr decltype(PentagonalPyramidal::size) PentagonalPyramidal::size;
constexpr decltype(PentagonalPyramidal::coordinates) PentagonalPyramidal::coordinates;
constexpr decltype(PentagonalPyramidal::rotations) PentagonalPyramidal::rotations;
constexpr decltype(PentagonalPyramidal::tetrahedra) PentagonalPyramidal::tetrahedra;


struct PentagonalBiPyramidal {
    /* 
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
     */
  static constexpr Symmetry::Name name = Symmetry::Name::PentagonalBiPyramidal;
  static constexpr unsigned size = 7;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
    if(a == b) {
      return 0;
    }

    if(a + b == 11) {
      return M_PI; // trans 5,6
    }

    if(ConstexprMagic::Math::XOR(a > 4, b > 4)) {
      return M_PI / 2; // any angle to axial index
    }

    // remainder are equatorial angles, like PentagonalPlanar
    unsigned absDiff = std::min(a - b, b - a);
    return std::min(
      absDiff,
      std::min(absDiff - 5, 5 - absDiff)
    ) * ConstexprMagic::Math::toRadians<double>(72);
  }
  static constexpr std::array<ConstexprMagic::Vector, 7> coordinates {{
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
    {4, 0, 1, 2, 3, 5, 6}, // C5 axial
    {1, 0, 4, 3, 2, 6, 5} // C2 equatorial on 3
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
    {{0, 1, 5, replaceMe}},
    {{1, 2, 5, replaceMe}},
    {{2, 3, 5, replaceMe}},
    {{3, 4, 5, replaceMe}},
    {{4, 0, 5, replaceMe}},
    {{0, 1, replaceMe, 6}},
    {{1, 2, replaceMe, 6}},
    {{2, 3, replaceMe, 6}},
    {{3, 4, replaceMe, 6}},
    {{4, 0, replaceMe, 6}}
  }};
#endif
};

constexpr decltype(PentagonalBiPyramidal::name) PentagonalBiPyramidal::name;
constexpr decltype(PentagonalBiPyramidal::size) PentagonalBiPyramidal::size;
constexpr decltype(PentagonalBiPyramidal::coordinates) PentagonalBiPyramidal::coordinates;
constexpr decltype(PentagonalBiPyramidal::rotations) PentagonalBiPyramidal::rotations;
constexpr decltype(PentagonalBiPyramidal::tetrahedra) PentagonalBiPyramidal::tetrahedra;


struct SquareAntiPrismatic {
  /*  Two representations, one oblique, the other more helpful.
   *  The first is a side-on view. 0 is mostly hidden by 1, 3 mostly hidden 
   *  by 2. 4 and 6 are in the viewing plane while 5 juts out above plane and
   *  7 dips behind plane. 
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
   * Angles:
   *
   *   in-plane cis (4, 5) -> 55°
   *   in-plane trans (4, 6) -> 148°
   *   short between planes (0, 4) -> 51°
   *   long between planes (0, 5) -> 175°
   *
   * from [ReF8]2-, Koz'min, P.A., Zhurnal Strukturnoi Khimii 1964 
   * HINT: use the ICSD for angle calculations and preview
   *
   */
  static constexpr Symmetry::Name name = Symmetry::Name::SquareAntiPrismatic;
  static constexpr unsigned size = 8;
  static const std::string stringName;
  static constexpr double angleFunction(const unsigned& a, const unsigned& b) {
#ifdef USE_CONSTEXPR_SQUARE_ANTIPRISMATIC_LOOKUP_TABLE
    if(a == b) {
      return 0;
    }

    return squareAntiprismaticAngles.at(
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
        return ConstexprMagic::Math::toRadians<double>(72.9875); 
      } 
      
      // otherwise trans
      return ConstexprMagic::Math::toRadians<double>(114.475);
    }
    
    // remaining cases are between planes
    unsigned minDiff = std::min(a - b, b - a);
    if(minDiff == 3 || minDiff == 4 || minDiff == 7) { // short
      return ConstexprMagic::Math::toRadians<double>(78.05);
    } 

    // last case is long between planes
    return ConstexprMagic::Math::toRadians<double>(142.275);
#endif
  }
  static constexpr std::array<ConstexprMagic::Vector, 8> coordinates {{
    {-0.00928803, 0.611568, 0.791137},
    {0.795627, 0.605641, -0.0132684},
    {0.795627, -0.605641, -0.0132684},
    {-0.00928803, -0.611568, 0.791137},
    {-0.396172, 0.852169, -0.341841},
    {0.293758, 0, -0.95588},
    {-0.396172, -0.852169, -0.341841},
    {-0.983087, 0, 0.183141}
  }};
  static constexpr std::array<
    std::array<unsigned, 8>,
    2
  > rotations {{
    {3, 0, 1, 2, 7, 4, 5, 6}, // C4 axial
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
    {5, 4, 7, 6, 1, 0, 3, 2}, 
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
    {{7, 0, 4, replaceMe}},
    {{0, 4, replaceMe, 1}},
    {{4, 1, 5, replaceMe}},
    {{1, 5, replaceMe, 2}},
    {{5, 2, 6, replaceMe}},
    {{2, 6, replaceMe, 3}},
    {{6, 3, 7, replaceMe}},
    {{3, 7, replaceMe, 0}}
  }};
#endif
};

constexpr decltype(SquareAntiPrismatic::name) SquareAntiPrismatic::name;
constexpr decltype(SquareAntiPrismatic::size) SquareAntiPrismatic::size;
constexpr decltype(SquareAntiPrismatic::coordinates) SquareAntiPrismatic::coordinates;
constexpr decltype(SquareAntiPrismatic::rotations) SquareAntiPrismatic::rotations;
constexpr decltype(SquareAntiPrismatic::tetrahedra) SquareAntiPrismatic::tetrahedra;


const std::string Linear::stringName {"linear"};
const std::string Bent::stringName {"bent"};
const std::string TrigonalPlanar::stringName {"trigonal planar"};
const std::string TrigonalPyramidal::stringName {"trigonal pyramidal"};
const std::string TShaped::stringName {"T-shaped"};
const std::string Tetrahedral::stringName {"tetrahedral"};
const std::string SquarePlanar::stringName {"square planar"};
const std::string Seesaw::stringName {"seesaw"};
const std::string TrigonalBiPyramidal::stringName {"trigonal bipyramidal"};
const std::string SquarePyramidal::stringName {"square pyramidal"};
const std::string PentagonalPlanar::stringName {"pentagonal planar"};
const std::string Octahedral::stringName {"octahedral"};
const std::string TrigonalPrismatic::stringName {"trigonal prismatic"};
const std::string PentagonalPyramidal::stringName {"pentagonal pyramidal"};
const std::string PentagonalBiPyramidal::stringName {"pentagonal bipyramidal"};
const std::string SquareAntiPrismatic::stringName {"square antiprismatic"};

//! Type collecting all types of the Symmetry classes.
using allSymmetryDataTypes = std::tuple<
  Linear, // 2
  Bent,
  TrigonalPlanar, // 3
  TrigonalPyramidal,
  TShaped,
  Tetrahedral, // 4
  SquarePlanar,
  Seesaw,
  SquarePyramidal, // 5
  TrigonalBiPyramidal,
  PentagonalPlanar,
  Octahedral, // 6
  TrigonalPrismatic,
  PentagonalPyramidal,
  PentagonalBiPyramidal, // 7
  SquareAntiPrismatic // 8
>;

/*! Template function implementation that helps to unpack a tuple of types and
 * forward them as template parameters to a template class implementing an "op()"
 * member. It groups the index sequence passed to it into a function parameter
 * pack and uses that pack in an expansion to extract all types contained in
 * the tuple, forwarding it to the pseudo functor.
 */
template<
  typename Tuple,
  template<typename ...> class PseudoFunctor,
  std::size_t... I
> constexpr auto unpackHelper(std::index_sequence<I...>) {
  return PseudoFunctor<
    std::tuple_element_t<
      I,
      Tuple
    >...
  >::op();
}

/*! Template function that unpacks a tuple of types and forwards them as
 * template parameters to a template class implementing an "op()" member.
 */
template<
  typename Tuple,
  template<typename ...> class PseudoFunctor
> constexpr auto unpackTupleToTemplateFunctor() {
  return unpackHelper<Tuple, PseudoFunctor>(
    std::make_index_sequence<
      std::tuple_size<Tuple>::value
    >()
  );
}

using AngleFunctionPtr = double(*)(const unsigned&, const unsigned&);

/*! Constructs an array of function pointers to all static angle functions
 * for runtime lookup
 */
template<typename ...SymmetryClasses>
struct angleFunctionFunctor {
  static constexpr std::array<AngleFunctionPtr, sizeof...(SymmetryClasses)> op() {
    std::array<AngleFunctionPtr, sizeof...(SymmetryClasses)> sizes = {{
      &SymmetryClasses::angleFunction...
    }};

    return sizes;
  }
};

//! Stub to find out the minimum angle returned in a specific symmetry class type
template<typename SymmetryClass> 
constexpr double smallestAngle() {
  double smallestAngle = M_PI;

  for(unsigned i = 0; i < SymmetryClass::size; ++i) {
    for(unsigned j = i + 1; j < SymmetryClass::size; ++j) {
      double returnedAngle = SymmetryClass::angleFunction(i, j);
      if(returnedAngle < smallestAngle) {
        smallestAngle = returnedAngle;
      }
    }
  }

  return smallestAngle;
}

/*! Functor to find out the minimum angle among all the symmetry class types
 * passed as template arguments
 */
template<typename ...SymmetryClasses>
struct minAngleFunctor {
  static constexpr double op() {
    const std::array<double, sizeof...(SymmetryClasses)> smallestAngles {{
      smallestAngle<SymmetryClasses>()...
    }};

    // C++17 min_element (isn't constexpr before)
    double minElement = smallestAngles.at(0);

    for(unsigned i = 1; i < sizeof...(SymmetryClasses); ++i) {
      if(smallestAngles.at(i) < minElement) {
        minElement = smallestAngles.at(i);
      }
    }

    return minElement;
  }
};

/*! Conversion function to make the dynamic rotations list type from the
 * constexpr data types given in a specifc symmetry class type
 */
template<size_t symmetrySize, size_t nRotations>
RotationsList makeRotations(
  const std::array<
    std::array<unsigned, symmetrySize>,
    nRotations
  >& constexprRotations
) {
  RotationsList rotations;

  for(const auto& rotation : constexprRotations) {
    rotations.emplace_back(
      rotation.begin(),
      rotation.end()
    );
  }

  return rotations;
}

/*! Conversion function to make the dynamic tetrahedron list type from the
 * constexpr data types given in a specifc symmetry class type
 */
template<size_t nTetrahedra>
TetrahedronList makeTetrahedra(
  const std::array<
    std::array<unsigned, 4>,
    nTetrahedra
  >& constexprTetrahedra
) {
  TetrahedronList tetrahedra;

  for(const auto& tetrahedron : constexprTetrahedra) {
    tetrahedra.push_back(
      TemplateMagic::map(
        tetrahedron,
        [](const unsigned& index) -> boost::optional<unsigned> {
          if(index == replaceMe) {
            return boost::none;
          }

          return index;
        }
      )
    );
  }

  return tetrahedra;
}

//! Conversion helper to Eigen type from constexpr vector type
Eigen::Vector3d toEigen(const ConstexprMagic::Vector& cVector) {
  return {
    cVector.data[0],
    cVector.data[1],
    cVector.data[2]
  };
}

/*! Conversion function to make the dynamic coordinates list type from the
 * constexpr data types given in a specifc symmetry class type
 */
template<size_t symmetrySize>
CoordinateList makeCoordinates(
  const std::array<ConstexprMagic::Vector, symmetrySize>& constexprCoordinates
) {
  return TemplateMagic::mapToVector(
    constexprCoordinates,
    toEigen
  );
}

/*! This constructs the SymmetryInformation instance from a specific symmetry
 * class type.
 */
template<typename SymmetryClass>
SymmetryInformation makeSymmetryInformation() {
  return {
    SymmetryClass::stringName,
    SymmetryClass::size,
    makeRotations(SymmetryClass::rotations),
    makeTetrahedra(SymmetryClass::tetrahedra),
    makeCoordinates(SymmetryClass::coordinates)
  };
}

/*! This creates a map initialization pair for a specific symmetry class type.
 * The key is the name, the mapped_type a SymmetryInformation instance
 */
template<typename SymmetryClass>
std::pair<Name, SymmetryInformation> makeMapInitPair() {
  return {
    SymmetryClass::name,
    makeSymmetryInformation<SymmetryClass>()
  };
}

/*! Creates the mapping between a symmetry class's name and it's dynamic
 * information in order to have runtime lookup based on symmetry names.
 */
template<typename ...SymmetryClasses>
struct symmetryInformationFunctor {
  static const std::map<Name, SymmetryInformation> op() {
    return {{
      makeMapInitPair<SymmetryClasses>()...
    }};
  };
};

constexpr auto angleFunctions = unpackTupleToTemplateFunctor<
  allSymmetryDataTypes,
  angleFunctionFunctor
>();

} // namespace data

const std::map<Name, SymmetryInformation> symmetryData = data::unpackTupleToTemplateFunctor<
  data::allSymmetryDataTypes,
  data::symmetryInformationFunctor
>();

// Shortcut functions
inline const std::string& name(const Name& name) {
  return symmetryData.at(name).stringName;
}

inline const unsigned& size(const Name& name) {
  return symmetryData.at(name).size;
}

inline const RotationsList& rotations(const Name& name) {
  return symmetryData.at(name).rotations;
}

inline data::AngleFunctionPtr angleFunction(const Name& name) {
  unsigned symmetryIndex = static_cast<unsigned>(name);
  return data::angleFunctions.at(symmetryIndex);
}

inline unsigned nameIndex(const Name& name) {
  return std::find(
    allNames.begin(),
    allNames.end(),
    name
  ) - allNames.begin();
}

inline const TetrahedronList& tetrahedra(const Name& name) {
  return symmetryData.at(name).tetrahedra;
}

// Derived data
constexpr double smallestAngle = data::unpackTupleToTemplateFunctor<
  data::allSymmetryDataTypes,
  data::minAngleFunctor
>();

} // namespace Symmetry

#endif
