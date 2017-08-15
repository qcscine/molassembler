#include "Properties.h"

#include <Eigen/Dense>

namespace Symmetry {

namespace properties {

std::vector<unsigned> applyRotation(
  const std::vector<unsigned>& indices,
  const Symmetry::Name& symmetryName,
  unsigned rotationFunctionIndex
) {
  std::vector<unsigned> retv;

  for(
    const auto& index : 
    Symmetry::rotations(symmetryName).at(rotationFunctionIndex)
  ) {
    retv.push_back(
      indices.at(index)
    );
  }
  
  return retv;
}

Eigen::Vector3d getCoordinates(
  const Symmetry::Name& symmetryName,
  const boost::optional<unsigned>& indexInSymmetryOption
) {
  assert(
    (
      indexInSymmetryOption 
      && indexInSymmetryOption.value() < Symmetry::size(symmetryName)
    ) || !indexInSymmetryOption
  );

  if(indexInSymmetryOption) {
    return Symmetry::symmetryData().at(symmetryName).coordinates.at(
      indexInSymmetryOption.value()
    );
  }

  return {0, 0, 0};
}

double getTetrahedronVolume(
  const Eigen::Vector3d& i,
  const Eigen::Vector3d& j,
  const Eigen::Vector3d& k,
  const Eigen::Vector3d& l
) {
  return (
    i - l
  ).dot(
    (j - l).cross(k - l)
  );
}

double calculateAngleDistortion(
  const Symmetry::Name& from,
  const Symmetry::Name& to,
  const std::vector<unsigned>& indexMapping
) {
  assert(indexMapping.size() == Symmetry::size(to));

  double angleDistortion = 0;

  for(unsigned i = 0; i < Symmetry::size(from); ++i) {
    for(unsigned j = i + 1; j < Symmetry::size(from); ++j) {
      angleDistortion += std::fabs(
        Symmetry::angleFunction(from)(i, j)
        - Symmetry::angleFunction(to)(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return angleDistortion;
}

boost::optional<unsigned> propagateIndexOptionalThroughMapping(
  const boost::optional<unsigned>& indexOptional,
  const std::vector<unsigned>& indexMapping
) {
  if(indexOptional) {
    return indexMapping.at(indexOptional.value());
  }

  return boost::none;
}


double calculateChiralDistortion(
  const Symmetry::Name& from,
  const Symmetry::Name& to,
  const std::vector<unsigned>& indexMapping
) {
  assert(indexMapping.size() == Symmetry::size(to));

  double chiralDistortion = 0;

  for(const auto& tetrahedron : Symmetry::tetrahedra(from)) {
    chiralDistortion += std::fabs(
      getTetrahedronVolume(
        getCoordinates(from, tetrahedron.at(0)),
        getCoordinates(from, tetrahedron.at(1)),
        getCoordinates(from, tetrahedron.at(2)),
        getCoordinates(from, tetrahedron.at(3))
      ) - getTetrahedronVolume(
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(0), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(1), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(2), indexMapping)
        ),
        getCoordinates(
          to,
          propagateIndexOptionalThroughMapping(tetrahedron.at(3), indexMapping)
        )
      )
    );
  }

  return chiralDistortion;
}


std::set<
  std::vector<unsigned>
> generateAllRotations(
  const Symmetry::Name& symmetryName,
  const std::vector<unsigned>& indices
) {
  using IndicesList = std::vector<unsigned>;

  std::set<IndicesList> allRotations = {indices};

  unsigned linkLimit = Symmetry::rotations(symmetryName).size();

  std::vector<unsigned> chain = {0};
  std::vector<IndicesList> chainStructures = {indices};

  // begin loop
  while(chain.front() < linkLimit) {
    // perform rotation
    // copy the last element in chainStructures
    auto generated = applyRotation(
      chainStructures.back(),
      symmetryName,
      chain.back()
    );

    // is it something new?
    if(allRotations.count(generated) == 0) {
      // add it to the set
      allRotations.insert(generated);

      // add it to the chain
      chainStructures.push_back(generated);
      chain.emplace_back(0);
    } else {
      // collapse the chain until we are at an incrementable position (if need be)
      while(
        chain.size() > 1
        && chain.back() == linkLimit - 1
      ) {
        chain.pop_back();
        chainStructures.pop_back();
      }

      ++chain.back();
    }
  }

  return allRotations;
}

std::vector<unsigned> symPosMapping(
  const std::vector<unsigned>& mapping
) {
  std::vector<unsigned> symmetryPositions (mapping.size());

  for(unsigned i = 0; i < mapping.size(); ++i) {
    symmetryPositions.at(
      mapping.at(i)
    ) = i;
  }
  
  return symmetryPositions;
}

DistortionInfo::DistortionInfo(
  const std::vector<unsigned>& passIndexMapping,
  const double& passTotalDistortion,
  const double& passChiralDistortion
) : indexMapping(passIndexMapping),
    totalDistortion(passTotalDistortion),
    chiralDistortion(passChiralDistortion)
{}

std::vector<DistortionInfo> ligandGainDistortions(
  const Symmetry::Name& symmetryFrom,
  const Symmetry::Name& symmetryTo
) {
  assert(Symmetry::size(symmetryTo) == Symmetry::size(symmetryFrom) + 1);

  std::vector<DistortionInfo> distortions;

  // Create a vector of indices for the new symmetry to use as mapping
  std::vector<unsigned> indexMapping (Symmetry::size(symmetryTo));
  std::iota(
    indexMapping.begin(),
    indexMapping.end(),
    0
  );

  std::set<
    std::vector<unsigned>
  > encounteredSymmetryMappings;

  do {
    if(encounteredSymmetryMappings.count(symPosMapping(indexMapping)) == 0) {
      distortions.emplace_back(
        indexMapping,
        calculateAngleDistortion(symmetryFrom, symmetryTo, indexMapping),
        calculateChiralDistortion(symmetryFrom, symmetryTo, indexMapping)
      );

      auto allRotations = generateAllRotations(
        symmetryTo,
        symPosMapping(indexMapping)
      );

      encounteredSymmetryMappings.insert(
        allRotations.begin(),
        allRotations.end()
      );
    }
  } while (std::next_permutation(indexMapping.begin(), indexMapping.end()));

  return distortions;
}

} // namespace properties

} // namespace Symmetry
