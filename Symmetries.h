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

#include "template_magic/TemplateMagic.h"

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

struct SymmetryInformation {
  const std::string stringName;
  const unsigned size;
  const RotationsList rotations;
  const AngleFunctionType angleFunction;
  const TetrahedronList tetrahedra;
  const CoordinateList coordinates;

  // Direct initialization
  SymmetryInformation(
    std::string&& stringName,
    unsigned&& size,
    RotationsList&& rotations,
    AngleFunctionType&& angleFunction,
    TetrahedronList&& tetrahedra,
    CoordinateList&& coordinates
  ) : stringName(stringName),
      size(size),
      rotations(rotations),
      angleFunction(angleFunction),
      tetrahedra(tetrahedra),
      coordinates(coordinates)
  {}
};

// Symmetry names list
enum class Name {
  Linear, // 2
  Bent,
  TrigonalPlanar, // 3
  TrigonalPyramidal,
  TShaped,
  Tetrahedral, // 4
  SquarePlanar,
  Seesaw,
  TrigonalBiPyramidal, // 5
  SquarePyramidal, 
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

extern const std::map<Name, SymmetryInformation> symmetryData;

namespace AngleFunctions {

using AngleFunctionPtr = double(*)(const unsigned&, const unsigned&);

struct SymmetrySizeAndAngles {
  constexpr SymmetrySizeAndAngles(
    const unsigned& size,
    double(*angleFunction)(const unsigned&, const unsigned&)
  ) : size(size),
      angleFunction(angleFunction) 
  {}

  const unsigned size;
  const AngleFunctionPtr angleFunction;
};

constexpr double linear(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  return M_PI;
}

constexpr double bent(const unsigned& a, const unsigned& b) {
  /* subject to a lot of variation, between 90 and 109 degrees pursuant to 
   * english wikipedia, using experimental data here to improve instances
   * of this geometry on e.g. O center would a big improvement to DG runs
   */
  if(a == b) {
    return 0;
  }

  return ConstexprMagic::Math::toRadians(107); 
}

constexpr double trigonalPlanar(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  return ConstexprMagic::Math::toRadians(120);
}

constexpr double trigonalPyramidal(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  return ConstexprMagic::Math::toRadians(107.5);
}

constexpr double tShaped(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  if((a + b) % 2 == 1) {
    return M_PI / 2;
  } 

  return M_PI;
}

constexpr double tetrahedral(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  return ConstexprMagic::Math::toRadians(109.5);
}

constexpr double squarePlanar(const unsigned& a, const unsigned& b) {
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

constexpr double seesaw(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  const auto& smaller = std::min(a, b);
  const auto& larger = std::max(a, b);
  if(smaller == 0 && larger == 3) {
    return M_PI;
  }

  if(smaller == 1 && larger == 2) {
    return ConstexprMagic::Math::toRadians(120);
  }

  return M_PI / 2;
}

constexpr double squarePyramidal(const unsigned& a, const unsigned& b) {
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

constexpr double trigonalBiPyramidal(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  unsigned smaller = std::min(a, b), larger = std::max(a, b);
  if(larger < 3) {
    // -> smaller < 2, this means either 0,1 0,2 1,2 axial
    return ConstexprMagic::Math::toRadians(120);
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

constexpr double pentagonalPlanar(const unsigned& a, const unsigned& b) {
  unsigned absDiff = std::min(a - b, b - a);
  return std::min(
    absDiff,
    std::min(absDiff - 5, 5 - absDiff)
  ) * ConstexprMagic::Math::toRadians(72);
}

constexpr double octahedral(const unsigned& a, const unsigned& b) {
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

constexpr double trigonalPrismatic(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  } 
  
  // Between plane symmetric
  if(std::min(a - b, b - a) == 3) {
    return ConstexprMagic::Math::toRadians(76);
  } 

  // In plane triangle
  if(
    (a < 3 && b < 3)
    || (a >= 3 && b >= 3)
  ) {
    return ConstexprMagic::Math::toRadians(86);
  } 

  // Between plane asymmetric
  return ConstexprMagic::Math::toRadians(134);
}

constexpr double pentagonalPyramidal(const unsigned& a, const unsigned& b) {
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
  ) * ConstexprMagic::Math::toRadians(72);
}

constexpr double pentagonalBiPyramidal(const unsigned& a, const unsigned& b) {
  if(a == b) {
    return 0;
  }

  if(a + b == 11) {
    return M_PI; // trans 5,6
  }

  if(TemplateMagic::XOR(a > 4, b > 4)) {
    return M_PI / 2; // any angle to axial index
  }

  // remainder are equatorial angles, like PentagonalPlanar
  unsigned absDiff = std::min(a - b, b - a);
  return std::min(
    absDiff,
    std::min(absDiff - 5, 5 - absDiff)
  ) * ConstexprMagic::Math::toRadians(72);
}

constexpr double squareAntiprismatic(const unsigned& a, const unsigned& b) {
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
      return ConstexprMagic::Math::toRadians(72.9875); 
    } 
    
    // otherwise trans
    return ConstexprMagic::Math::toRadians(114.475);
  }
  
  // remaining cases are between planes
  unsigned minDiff = std::min(a - b, b - a);
  if(minDiff == 3 || minDiff == 4 || minDiff == 7) { // short
    return ConstexprMagic::Math::toRadians(78.05);
  } 

  // last case is long between planes
  return ConstexprMagic::Math::toRadians(142.275);
#endif
}

constexpr std::array<SymmetrySizeAndAngles, nSymmetries> sizeAndAngles {
SymmetrySizeAndAngles {2, &AngleFunctions::linear},
  SymmetrySizeAndAngles {2, &AngleFunctions::bent},
  SymmetrySizeAndAngles {3, &AngleFunctions::trigonalPlanar},
  SymmetrySizeAndAngles {3, &AngleFunctions::trigonalPyramidal},
  SymmetrySizeAndAngles {3, &AngleFunctions::tShaped},
  SymmetrySizeAndAngles {4, &AngleFunctions::tetrahedral},
  SymmetrySizeAndAngles {4, &AngleFunctions::squarePlanar},
  SymmetrySizeAndAngles {4, &AngleFunctions::seesaw},
  SymmetrySizeAndAngles {5, &AngleFunctions::squarePyramidal},
  SymmetrySizeAndAngles {5, &AngleFunctions::trigonalBiPyramidal},
  SymmetrySizeAndAngles {5, &AngleFunctions::pentagonalPlanar},
  SymmetrySizeAndAngles {6, &AngleFunctions::octahedral},
  SymmetrySizeAndAngles {6, &AngleFunctions::trigonalPrismatic},
  SymmetrySizeAndAngles {6, &AngleFunctions::pentagonalPyramidal},
  SymmetrySizeAndAngles {7, &AngleFunctions::pentagonalBiPyramidal},
  SymmetrySizeAndAngles {8, &AngleFunctions::squareAntiprismatic}
};

constexpr unsigned getIndexOfName(const Name& name) {
  unsigned i = 0;

  while(allNames.at(i) != name) {
    ++i;
  }

  return i;
}

constexpr double minAngle() {
  double smallestAngle = M_PI;

  for(unsigned i = 0; i < nSymmetries; i++) {
    auto& name = allNames.at(i);
    unsigned symmetryIndex = getIndexOfName(name);
    unsigned maxIndex = sizeAndAngles.at(symmetryIndex).size;

    for(unsigned i = 0; i < maxIndex; ++i) {
      for(unsigned j = i + 1; j < maxIndex; ++j) {
        double returnedAngle = sizeAndAngles.at(symmetryIndex).angleFunction(i, j);

        if(returnedAngle < smallestAngle) {
          smallestAngle = returnedAngle;
        }
      }
    }
  }

  return smallestAngle;
}

} // namespace AngleFunctions

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

inline const AngleFunctionType& angleFunction(const Name& name) {
  return symmetryData.at(name).angleFunction;
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
constexpr double smallestAngle = AngleFunctions::minAngle();

} // namespace Symmetry

#endif
