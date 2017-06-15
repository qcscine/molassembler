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
 * - C++17 improvements:
 *   - Make angle calculation constexpr so that smallestAngle is also constexpr
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
extern const std::vector<Name> allNames;
extern const std::map<Name, SymmetryInformation> symmetryData;

// derived data
extern const double smallestAngle;

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

} // namespace Symmetry

#endif
