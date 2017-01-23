#ifndef LIB_INCLUDE_SYMMETRIES_H
#define LIB_INCLUDE_SYMMETRIES_H

#include <map>
#include <vector>
#include <functional>

namespace Symmetry {

enum class Name {
  Linear,
  TrigonalPlanar,
  Tetrahedral,
  SquarePlanar,
  SquarePyramidal,
  TrigonalBiPyramidal,
  Octahedral
};

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

extern const std::map<Name, TupleType> symmetryData;

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
