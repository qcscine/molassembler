#include "Symmetries.h"

#include "constexpr_magic/Containers.h"
#include "constexpr_magic/DynamicSet.h"
#include "constexpr_magic/Array.h"

/*! @file
 *
 * Constexpr parallel to Properties.h. Contains a slew of computations on the
 * base symmetry data to compute derived properties. You can e.g. apply
 * rotations, extract rotation multiplicities, calculate angular and chiral
 * distortions between pairs of symmetries, etc.
 */

/* TODO
 * - What is the exact difference between propagateIndexMapping and
 *   symPosMapping?
 */

namespace Symmetry {

/* Typedef to use ConstexprMagic's Array instead of std::array as the underlying
 * base array type since C++14's std::array has too few members marked constexpr
 * as to be useful. When C++17 rolls around, replace this with std::array!
 */
template<typename T, size_t size> 
using ArrayType = ConstexprMagic::Array<T, size>;

namespace detail {

// Iota helper, required to expand the index sequence generated in base
template<size_t size, size_t... Indices>
constexpr ArrayType<unsigned, size> iotaImpl(
  std::index_sequence<Indices...>
) {
  return {
    Indices...
  };
}

} // namespace detail

template<size_t size>
constexpr ArrayType<unsigned, size> iota() {
  return detail::iotaImpl<size>(
    std::make_index_sequence<size>{}
  );
}

// applyRotation helper function
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

//! Applies a symmetry group rotation to an array of indices
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
  } 

  return rotationMultiplicityImpl<SymmetryClass, rotationFunctionIndex>(
    applyRotation<SymmetryClass>(
      runningIndices,
      rotationFunctionIndex
    ),
    count + 1
  );
}

/*!
 * Calculates the multiplicity of a symmetry group's rotation specified via an
 * index in that symmetry's list of rotations
 */
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

//! Calculates all multiplicities of a symmetry group's rotations.
template<typename SymmetryClass>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()> rotationMultiplicities() {
  return rotationMultiplicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
}

/*!
 * Template function generating all rotation multiplicities for a specified
 * symmetry group type
 */
template<typename SymmetryClass>
struct allRotationMultiplicities {
  static constexpr ArrayType<
    unsigned,
    SymmetryClass::rotations
  > value = rotationMultiplicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
};

/*!
 * Fetches the coordinates of an index in a Symmetry, properly handling the
 * boost::none -> origin mapping.
 */
template<typename SymmetryClass>
constexpr ConstexprMagic::Vector getCoordinates(const unsigned& indexInSymmetry) {
  if(indexInSymmetry != ORIGIN_PLACEHOLDER) {
    return SymmetryClass::coordinates.at(indexInSymmetry);
  }

  return {0, 0, 0};
}

/*!
 * Calculates the volume of a tetrahedron spanned by four positions in
 * three-dimensional space.
 */
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

/*!
 * Calculates the angular distortion between a pair of symmetries and a given
 * index mapping between them that specifies how symmetry positions are mapped
 * between the two.
 */
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
  double angularDistortion = 0;

  for(unsigned i = 0; i < SymmetryClassFrom::size; ++i) {
    for(unsigned j = i + 1; j < SymmetryClassFrom::size; ++j) {
      angularDistortion += ConstexprMagic::Math::abs(
        SymmetryClassFrom::angleFunction(i, j)
        - SymmetryClassTo::angleFunction(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return angularDistortion;
}

/*!
 * Propagates a source symmetry position through an index mapping, properly
 * handling the origin placeholder special unsigned value.
 */
template<size_t size>
constexpr unsigned propagateSymmetryPosition(
  const unsigned& symmetryPosition,
  const ArrayType<unsigned, size>& indexMapping
) {
  if(symmetryPosition != ORIGIN_PLACEHOLDER) {
    return indexMapping.at(symmetryPosition);
  }

  return ORIGIN_PLACEHOLDER;
}

/*!
 * Calculate the chiral distortion between the source symmetry and the target
 * symmetry specified by an index mapping.
 */
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

/*!
 * Calculate the maximum number of non-superimposable rotations a symmetry
 * class can produce for entirely unequal indices
 */
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

// Some helper types
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
  maxRotations<SymmetryClass>() * 2 // TODO factor is entirely arbitrary
>;

template<typename SymmetryClass>
using ChainArrayType = ConstexprMagic::DynamicArray<
  unsigned,
  maxRotations<SymmetryClass>() * 2 // TODO factor is entirely arbitrary
>;

//! Generates all rotations of a sequence of indices within a symmetry group
template<typename SymmetryClass>
constexpr auto generateAllRotations(const IndicesList<SymmetryClass>& indices) {
  RotationsSetType<SymmetryClass> rotations;
  rotations.insert(indices);

  ChainStructuresArrayType<SymmetryClass> chainStructures;
  chainStructures.push_back(indices);

  ChainArrayType<SymmetryClass> chain {0u};

  // The very last rotation isn't found for PentagonalPyramidal for some reason
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
        chain.size() > 1 
        && chain.back() == SymmetryClass::rotations.size() - 1
      ) {
        chain.pop_back();
        chainStructures.pop_back();
      }

      // increment
      ++chain.back();
    }
  }

  return rotations;
}

/*! 
 * Data struct to collect the results of calculating the ideal index mappings
 * between pairs of indices
 */
template<typename SymmetryClass>
struct MappingsReturnType {
  using MappingsList = ConstexprMagic::DynamicSet<
    ArrayType<unsigned, SymmetryClass::size>,
    20
  >;

  MappingsList mappings;
  double angularDistortion, chiralDistortion;

  constexpr MappingsReturnType(
    MappingsList&& mappings,
    double&& angularDistortion,
    double&& chiralDistortion
  ) : mappings(mappings),
      angularDistortion(angularDistortion),
      chiralDistortion(chiralDistortion)
  {}
};

/*!
 * Calculates the ideal index mappings when adding a ligand to a particular
 * symmetry.
 *
 * TODO generalize to more equal size target symmetry and perhaps also ligand
 * loss as soon as dynamic variant works
 */
template<class SymmetryClassFrom, class SymmetryClassTo>
constexpr auto ligandGainMappings() {
  static_assert(
    SymmetryClassTo::size == SymmetryClassFrom::size + 1,
    "Ligand gain pathway calculation must involve symmetry size increase"
  );

  using IndexMappingType = ArrayType<unsigned, SymmetryClassTo::size>;

  ConstexprMagic::DynamicSet<IndexMappingType, 20> bestMappings;
  double lowestAngleDistortion = 100;
  double lowestChiralDistortion = 100;

  auto indexMapping = iota<SymmetryClassTo::size>();

  ConstexprMagic::DynamicSet<
    IndexMappingType,
    ConstexprMagic::Math::factorial(SymmetryClassTo::size)
  > encounteredMappings;

  do {
    if(!encounteredMappings.contains(symPosMapping(indexMapping))) {
      auto angularDistortion = calculateAngleDistortion<
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
       * The boolean cases below are, AFAICT, the fewest comparisons.
       * Feel free to see if you can find some way with fewer.
       */

      bool addMapping = (
        angularDistortion < lowestAngleDistortion
        || (
          angularDistortion == lowestAngleDistortion
          && chiralDistortion <= lowestChiralDistortion
        )
      );

      bool clearExisting = (
        addMapping
        && !(
          angularDistortion == lowestAngleDistortion
          && chiralDistortion == lowestChiralDistortion
        )
      );
        
      if(clearExisting) {
        bestMappings.clear();
        lowestAngleDistortion = angularDistortion;
        lowestChiralDistortion = chiralDistortion;
      }

      if(addMapping) {
        bestMappings.insert(indexMapping);
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
