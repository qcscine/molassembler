#ifndef LIB_INCLUDE_SYMMETRIES_H
#define LIB_INCLUDE_SYMMETRIES_H

#include <map>
#include <vector>
#include <functional>
#include <algorithm>
#include "boost/optional.hpp"
#include <Eigen/Core>

/* TODO
 * - Debug and Release builds
 * - Square antiprismatic coordinates / angles need improvement, see tests.
 * - Improve trigonal pyramidal to get 107.5 angle as a parameter. Currently, 
 *   the rotation angle choice of 111.5 works well, but not ideal.
 * - Consider making constexpr calculation of angles from coordinates into
 *   const lookup table (coordinates would need to be std::array<double, 3>),
 *   and would require a constexpr math library
 */

namespace Symmetry {

/* Typedefs */
using RotationsType = std::vector<
  std::vector<unsigned>
>;

/* All angle functions can be called with arbitrary (valid) parameters
 * without failing. Valid here means that a != b and less than the size of
 * the symmetry requested.
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
  const RotationsType rotations;
  const AngleFunctionType angleFunction;
  const TetrahedronList tetrahedra;
  const CoordinateList coordinates;

  // Direct initialization
  SymmetryInformation(
    std::string&& stringName,
    unsigned&& size,
    RotationsType&& rotations,
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

// Shortcut functions
inline const std::string& name(const Name& name) {
  return symmetryData.at(name).stringName;
}

inline const unsigned& size(const Name& name) {
  return symmetryData.at(name).size;
}

inline const RotationsType& rotations(const Name& name) {
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

} // eo namespace

#endif
