#ifndef LIB_INCLUDE_SYMMETRIES_H
#define LIB_INCLUDE_SYMMETRIES_H

#include <map>
#include <vector>
#include <functional>

/* TODO
 */

namespace Symmetry {

/* Typedefs */
using RotationsType = std::vector<
  std::vector<unsigned>
>;

using AngleFunctionType = std::function<
  double(const unsigned&, const unsigned&)
>;

using TupleType = std::tuple<
  std::string,
  unsigned,
  RotationsType,
  AngleFunctionType
>;

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
extern const std::map<Name, TupleType> symmetryData;

// getter functions
const std::string& name(const Name& name) {
  return std::get<0>(
    symmetryData.at(name)
  );
}

unsigned size(const Name& name) {
  return std::get<1>(
    symmetryData.at(name)
  );
}

const RotationsType& rotations(const Name& name) {
  return std::get<2>(
    symmetryData.at(name)
  );
}

const AngleFunctionType& angleFunction(const Name& name) {
  return std::get<3>(
    symmetryData.at(name)
  );
}

} // eo namespace

#endif
