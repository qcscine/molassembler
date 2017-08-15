#ifndef INCLUDE_SYMMETRY_DYNAMIC_PROPERTIES_CALCULATION_H
#define INCLUDE_SYMMETRY_DYNAMIC_PROPERTIES_CALCULATION_H

#include <vector>
#include <numeric>

#include "Symmetries.h"

namespace Symmetry {

namespace properties {

std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const Symmetry::Name& symmetryName,
  unsigned rotationFunctionIndex
);

Eigen::Vector3d getCoordinates(
  const Symmetry::Name& symmetryName,
  const boost::optional<unsigned>& indexInSymmetryOption
);

double getTetrahedronVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
);

double calculateAngleDistortion(
  const Symmetry::Name& from,
  const Symmetry::Name& to,
  const std::vector<unsigned>& indexMapping
);

boost::optional<unsigned> propagateIndexOptionalThroughMapping(
  const boost::optional<unsigned>& indexOptional,
  const std::vector<unsigned>& indexMapping
);


double calculateChiralDistortion(
  const Symmetry::Name& from,
  const Symmetry::Name& to,
  const std::vector<unsigned>& indexMapping
);


std::set<
  std::vector<unsigned>
> generateAllRotations(
  const Symmetry::Name& symmetryName,
  const std::vector<unsigned>& indices
);

std::vector<unsigned> symPosMapping(const std::vector<unsigned>& mapping);

struct DistortionInfo {
  std::vector<unsigned> indexMapping;
  double totalDistortion;
  double chiralDistortion;

  DistortionInfo(
    const std::vector<unsigned>& passIndexMapping,
    const double& passTotalDistortion,
    const double& passChiralDistortion
  );
};

std::vector<DistortionInfo> ligandGainDistortions(
  const Symmetry::Name& symmetryFrom,
  const Symmetry::Name& symmetryTo
);

} // namespace properties

namespace detail {

template<typename NumericType>
std::vector<NumericType> iota(NumericType nValues) {
  std::vector<NumericType> values (nValues);

  std::iota(
    values.begin(),
    values.end(),
    NumericType {0}
  );

  return values;
}

} // namespace detail

} // namespace Symmetry

#endif
