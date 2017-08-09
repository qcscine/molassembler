#include "Symmetries.h"

#include "boost/hana.hpp"
#include "constexpr_magic/Containers.h"
#include "constexpr_magic/DynamicSet.h"
#include "constexpr_magic/Math.h"
#include "constexpr_magic/Array.h"

namespace Symmetry {

// When C++17 rolls around, replace this with std::array!
template<typename T, size_t size> 
using ArrayType = ConstexprMagic::Array<T, size>;

template<size_t size, size_t... Indices>
constexpr ArrayType<unsigned, size> iotaImpl(
  std::index_sequence<Indices...>
) {
  return {
    Indices...
  };
}

template<size_t size>
constexpr ArrayType<unsigned, size> iota() {
  return iotaImpl<size>(
    std::make_index_sequence<size>{}
  );
}

template<typename SymmetryClass, size_t... Indices>
constexpr ArrayType<unsigned, SymmetryClass::size> applyRotationImpl(
  const ArrayType<unsigned, SymmetryClass::size>& indices, 
  const unsigned& rotationFunctionIndex,
  std::index_sequence<Indices...>
) {
  return {
    indices.at(
      SymmetryClass::rotations.at(rotationFunctionIndex).at(Indices)
    )...
  };
}

template<typename SymmetryClass>
constexpr ArrayType<unsigned, SymmetryClass::size> applyRotation(
  const ArrayType<unsigned, SymmetryClass::size>& indices, 
  const unsigned& rotationFunctionIndex
) {
  return applyRotationImpl<SymmetryClass>(
    indices,
    rotationFunctionIndex,
    std::make_index_sequence<SymmetryClass::size>{}
  );
}

template<typename SymmetryClass>
constexpr unsigned rotationMultiplicityImpl(
  const unsigned& rotationFunctionIndex,
  const ArrayType<unsigned, SymmetryClass::size>& runningIndices,
  const unsigned& count
) {
  if(
    ConstexprMagic::arraysEqual(
      runningIndices,
      iota<SymmetryClass::size>()
    )
  ) {
    return count;
  } else {
    return rotationMultiplicityImpl<SymmetryClass>(
      rotationFunctionIndex,
      applyRotation<SymmetryClass>(
        runningIndices,
        rotationFunctionIndex
      ),
      count + 1
    );
  }
}

template<typename SymmetryClass>
constexpr unsigned rotationMultiplicity(
  const unsigned& rotationFunctionIndex
) {
  return rotationMultiplicityImpl<SymmetryClass>(
    rotationFunctionIndex,
    applyRotation<SymmetryClass>(
      iota<SymmetryClass::size>(),
      rotationFunctionIndex
    ),
    1
  );
}

template<typename SymmetryClass, unsigned rotationFunctionIndex>
constexpr unsigned rotationMultiplicityImpl(
  const ArrayType<unsigned, SymmetryClass::size>& runningIndices,
  const unsigned& count
) {
  if(
    ConstexprMagic::arraysEqual(
      runningIndices,
      iota<SymmetryClass::size>()
    )
  ) {
    return count;
  } else {
    return rotationMultiplicityImpl<SymmetryClass, rotationFunctionIndex>(
      applyRotation<SymmetryClass>(
        runningIndices,
        rotationFunctionIndex
      ),
      count + 1
    );
  }
}

template<typename SymmetryClass, unsigned rotationFunctionIndex>
constexpr unsigned rotationMultiplicity() {
  return rotationMultiplicityImpl<SymmetryClass, rotationFunctionIndex>(
    applyRotation<SymmetryClass>(
      iota<SymmetryClass::size>(),
      rotationFunctionIndex
    ),
    1
  );
}

template<typename SymmetryClass, size_t ...Inds>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()> 
rotationMultiplicitiesImpl(std::index_sequence<Inds...>) {
  return {
    rotationMultiplicity<SymmetryClass, Inds>()...
  };
}

template<typename SymmetryClass>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()> rotationMultiplicities() {
  return rotationMultiplicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
}

template<typename SymmetryClass>
struct allRotationMultiplicities {
  static constexpr ArrayType<
    unsigned,
    SymmetryClass::rotations
  > value = rotationMultiplicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
};


template<typename SymmetryClass>
constexpr ConstexprMagic::Vector getCoordinates(const unsigned& indexInSymmetry) {
  if(indexInSymmetry != replaceMe) {
    return SymmetryClass::coordinates.at(indexInSymmetry);
  }

  return {0, 0, 0};
}

constexpr double getTetrahedronVolume(
  const ConstexprMagic::Vector& i,
  const ConstexprMagic::Vector& j,
  const ConstexprMagic::Vector& k,
  const ConstexprMagic::Vector& l
) {
  return (i - l).dot(
    (j - l).cross(k - l)
  );
}

template<typename SymmetryClassFrom, typename SymmetryClassTo>
constexpr double calculateAngleDistortion(
  const ArrayType<
    unsigned,
    ConstexprMagic::Math::max(
      SymmetryClassFrom::size,
      SymmetryClassTo::size
    )
  >& indexMapping
) {
  double angleDistortion = 0;

  for(unsigned i = 0; i < SymmetryClassFrom::size; ++i) {
    for(unsigned j = i + 1; j < SymmetryClassFrom::size; ++j) {
      angleDistortion += ConstexprMagic::Math::abs(
        SymmetryClassFrom::angleFunction(i, j)
        - SymmetryClassTo::angleFunction(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return angleDistortion;
}

template<size_t size>
constexpr unsigned propagateSymmetryPosition(
  const unsigned& symmetryPosition,
  const ArrayType<unsigned, size>& indexMapping
) {
  if(symmetryPosition != replaceMe) {
    return indexMapping.at(symmetryPosition);
  }

  return replaceMe;
}

template<typename SymmetryClassFrom, typename SymmetryClassTo>
constexpr double calculateChiralDistortion(
  const ArrayType<
    unsigned,
    ConstexprMagic::Math::max(
      SymmetryClassFrom::size,
      SymmetryClassTo::size
    )
  >& indexMapping
) {
  double chiralDistortion = 0;

  // C++17:
  // for(const auto& tetrahedron : SymmetryClassFrom::tetrahedra) {
  for(unsigned i = 0; i < SymmetryClassFrom::tetrahedra.size(); ++i) {
    const auto& tetrahedron = SymmetryClassFrom::tetrahedra.at(i);

    chiralDistortion += ConstexprMagic::Math::abs(
      getTetrahedronVolume(
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(0)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(1)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(2)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(3))
      ) - getTetrahedronVolume(
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(0), indexMapping)),
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(1), indexMapping)),
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(2), indexMapping)),
        getCoordinates<SymmetryClassTo>(propagateSymmetryPosition(tetrahedron.at(3), indexMapping))
      )
    );
  }

  return chiralDistortion;
}

template<size_t size>
constexpr ArrayType<unsigned, size> symPosMapping(
  const ArrayType<unsigned, size>& mapping
) {
  ArrayType<unsigned, size> symmetryPositions;

  for(unsigned i = 0; i < size; ++i) {
    symmetryPositions.at(
      mapping.at(i)
    ) = i;
  }

  return symmetryPositions;
}

template<size_t size>
struct DistortionInfo {
  ArrayType<unsigned, size> indexMapping;
  double totalDistortion;
  double chiralDistortion;

  constexpr DistortionInfo(
    const ArrayType<unsigned, size>& passIndexMapping,
    const double& passTotalDistortion,
    const double& passChiralDistortion
  ) : indexMapping(passIndexMapping),
      totalDistortion(passTotalDistortion),
      chiralDistortion(passChiralDistortion)
  {}

  constexpr DistortionInfo()
    : indexMapping(iota<size>()),
      totalDistortion(100),
      chiralDistortion(100)
  {}
};


template<typename SymmetryClass>
constexpr unsigned maxRotations() {
  /* For each rotation in the symmetry class, figure out the multiplicity, i.e.
   * how often a rotation has to be applied to return to identity
   */
  constexpr auto symmetryRotationMultiplicities = rotationMultiplicities<SymmetryClass>();

  return ConstexprMagic::reduce(
    symmetryRotationMultiplicities,
    1u,
    std::multiplies<unsigned>()
  );
}


template<typename SymmetryClass>
using IndicesList = ArrayType<unsigned, SymmetryClass::size>;

template<typename SymmetryClass>
using RotationsSetType = ConstexprMagic::DynamicSet<
  IndicesList<SymmetryClass>, 
  maxRotations<SymmetryClass>()
>;

template<typename SymmetryClass>
using ChainStructuresArrayType = ConstexprMagic::DynamicArray<
  IndicesList<SymmetryClass>,
  SymmetryClass::rotations.size() * 3
>;

template<typename SymmetryClass>
using ChainArrayType = ConstexprMagic::DynamicArray<
  unsigned,
  SymmetryClass::rotations.size() * 3
>;


template<typename SymmetryClass>
constexpr void generateAllRotationsImpl(
  RotationsSetType<SymmetryClass>& rotations,
  ChainStructuresArrayType<SymmetryClass>& chainStructures,
  ChainArrayType<SymmetryClass>& chain
) {
  while(
    chain.front() < SymmetryClass::rotations.size() 
    && rotations.size() < maxRotations<SymmetryClass>()
  ) {

    auto generated = applyRotation<SymmetryClass>(
      chainStructures.back(),
      chain.back()
    );

    if(!rotations.contains(generated)) {
      rotations.insert(generated);
      chainStructures.push_back(generated);
      chain.push_back(0);
    } else {
      // collapse the chain until we are at an incrementable position (if need be)
      while(
        chain.size() > 0 
        && chain.back() == SymmetryClass::rotations.size() - 1
      ) {
        chain.pop_back();
        chainStructures.pop_back();
      }

      // increment
      ++chain.back();
    }
  }

}

template<typename SymmetryClass>
constexpr auto generateAllRotations(const IndicesList<SymmetryClass>& indices) {
  RotationsSetType<SymmetryClass> rotations;
  rotations.insert(indices);

  ChainStructuresArrayType<SymmetryClass> chainStructures;
  chainStructures.push_back(indices);

  ChainArrayType<SymmetryClass> chain {0u};

  generateAllRotationsImpl<SymmetryClass>(
    rotations,
    chainStructures,
    chain
  );

  return rotations;
}

template<typename SymmetryClass>
struct MappingsReturnType {
  using MappingsList = ConstexprMagic::DynamicArray<
    ArrayType<unsigned, SymmetryClass::size>,
    20
  >;

  MappingsList mappings;
  double angleDistortion, chiralDistortion;

  constexpr MappingsReturnType(
    MappingsList&& mappings,
    double&& angleDistortion,
    double&& chiralDistortion
  ) : mappings(mappings),
      angleDistortion(angleDistortion),
      chiralDistortion(chiralDistortion)
  {}
};

template<class SymmetryClassFrom, class SymmetryClassTo>
constexpr auto ligandGainMappings() {
  static_assert(
    SymmetryClassTo::size == SymmetryClassFrom::size + 1,
    "Ligand gain pathway calculation must involve symmetry size increase"
  );

  using IndexMappingType = ArrayType<unsigned, SymmetryClassTo::size>;

  ConstexprMagic::DynamicArray<IndexMappingType, 20> bestMappings;
  double lowestAngleDistortion = 100;
  double lowestChiralDistortion = 100;

  auto indexMapping = iota<SymmetryClassTo::size>();

  ConstexprMagic::DynamicSet<
    IndexMappingType,
    ConstexprMagic::Math::factorial(SymmetryClassTo::size)
  > encounteredMappings;

  do {
    if(!encounteredMappings.contains(symPosMapping(indexMapping))) {
      auto angleDistortion = calculateAngleDistortion<
        SymmetryClassFrom,
        SymmetryClassTo
      >(indexMapping);

      auto chiralDistortion = calculateChiralDistortion<
        SymmetryClassFrom,
        SymmetryClassTo
      >(indexMapping);

      /* Summary of cases:
       * - If any improvement is made on angular distortion, clear all and add
       * - If angular distortion is equal but chiral distortion is improved, 
       *   clear all and add
       * - If angular distortion and chiral distortion are equal, just add
       *
       * The boolean cases below are, AFAICT, the least amount of comparisons.
       * Feel free to see if you can find some way with fewer.
       */

      bool addMapping = (
        angleDistortion < lowestAngleDistortion
        || (
          angleDistortion == lowestAngleDistortion
          && chiralDistortion <= lowestChiralDistortion
        )
      );

      bool clearExisting = (
        addMapping
        && !(
          angleDistortion == lowestAngleDistortion
          && chiralDistortion == lowestChiralDistortion
        )
      );
        
      if(clearExisting) {
        bestMappings.clear();
        lowestAngleDistortion = angleDistortion;
        lowestChiralDistortion = chiralDistortion;
      }

      if(addMapping) {
        bestMappings.push_back(indexMapping);
      }

      // Add all rotations to the encountered mappings
      auto allRotations = generateAllRotations<SymmetryClassTo>(
        symPosMapping(indexMapping)
      );

      for(const auto& rotation : allRotations) {
        encounteredMappings.insert(rotation);
      }
    }
  } while(ConstexprMagic::inPlaceNextPermutation(indexMapping));

  return MappingsReturnType<SymmetryClassTo>(
    std::move(bestMappings),
    std::move(lowestAngleDistortion),
    std::move(lowestChiralDistortion)
  );
}

} // namespace Symmetry
