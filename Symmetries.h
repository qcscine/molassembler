#ifndef LIB_INCLUDE_SYMMETRIES_H
#define LIB_INCLUDE_SYMMETRIES_H

#include <map>
#include <vector>
#include <functional>
#include <algorithm>

/* TODO
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

struct SymmetryInformation {
  const std::string stringName;
  const unsigned size;
  const RotationsType rotations;
  const AngleFunctionType angleFunction;

  // Direct initialization
  SymmetryInformation(
    std::string&& stringName,
    unsigned&& size,
    RotationsType&& rotations,
    AngleFunctionType&& angleFunction
  ) : stringName(stringName),
      size(size),
      rotations(rotations),
      angleFunction(angleFunction)
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
inline const std::string& name(const Name& name);
inline unsigned size(const Name& name);
inline const RotationsType& rotations(const Name& name);
inline const AngleFunctionType& angleFunction(const Name& name);

inline const std::string& name(const Name& name) {
  return symmetryData.at(name).stringName;
}

inline unsigned size(const Name& name) {
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

// helpers
unsigned nameIndex(const Name& name);

} // eo namespace

#endif
