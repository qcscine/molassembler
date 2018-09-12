#ifndef INCLUDE_SYMMETRY_INFORMATION_CONSTEXPR_PROPERTIES_H
#define INCLUDE_SYMMETRY_INFORMATION_CONSTEXPR_PROPERTIES_H

#include "chemical_symmetries/Symmetries.h"

#include "temple/constexpr/Array.h"
#include "temple/constexpr/Bitset.h"
#include "temple/constexpr/DynamicSet.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Cache.h"

/*! @file
 *
 * @brief Compile-time symmetry property calculations
 *
 * Constexpr parallel to Properties.h. Contains a slew of computations on the
 * base symmetry data to compute derived properties. You can e.g. apply
 * rotations, extract rotation multiplicities, calculate angular and chiral
 * distortions between pairs of symmetries, etc.
 */

/* TODO
 * - remove hashIndexList in favor of the bitset hashing algorithm for
 *   constexpr speed
 * - remove unused stuff
 */

namespace Symmetry {

namespace constexprProperties {

constexpr double floatingPointEqualityTolerance = 1e-4;

//! Stub to find out the minimum angle returned in a specific symmetry class type
template<typename SymmetryClass>
constexpr double calculateSmallestAngle() {
  double smallestAngle = M_PI;

  for(unsigned i = 0; i < SymmetryClass::size; ++i) {
    for(unsigned j = i + 1; j < SymmetryClass::size; ++j) {
      double returnedAngle = SymmetryClass::angleFunction(i, j);
      if(returnedAngle < smallestAngle) {
        smallestAngle = returnedAngle;
      }
    }
  }

  return smallestAngle;
}

/*! Functor to find out the minimum angle among all the symmetry class types
 * passed as template arguments
 */
template<typename ...SymmetryClasses>
struct minAngleFunctor {
  static constexpr double value() {
    const std::array<double, sizeof...(SymmetryClasses)> smallestAngles {{
      calculateSmallestAngle<SymmetryClasses>()...
    }};

    // C++17 min_element (isn't constexpr before)
    double minElement = smallestAngles.at(0);

    for(unsigned i = 1; i < sizeof...(SymmetryClasses); ++i) {
      if(smallestAngles.at(i) < minElement) {
        minElement = smallestAngles.at(i);
      }
    }

    return minElement;
  }
};

/*!
 * Typedef to use temple's Array instead of std::array as the underlying
 * base array type since C++14's std::array has too few members marked constexpr
 * as to be useful. When C++17 rolls around, replace this with std::array!
 */
template<typename T, size_t size>
using ArrayType = temple::Array<T, size>;

//! Generate an integer sequence to use with stereopermutations
template<typename SymmetryClass>
constexpr auto startingIndexSequence() {
  return temple::iota<
    ArrayType,
    unsigned,
    SymmetryClass::size
  >();
}

//! Helper to perform applyRotation in constexpr fashion
template<typename SymmetryClass, size_t... Indices>
constexpr ArrayType<unsigned, SymmetryClass::size> applyRotationImpl(
  const ArrayType<unsigned, SymmetryClass::size>& indices,
  const unsigned rotationFunctionIndex,
  std::index_sequence<Indices...> /* inds */
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
  const unsigned rotationFunctionIndex
) {
  return applyRotationImpl<SymmetryClass>(
    indices,
    rotationFunctionIndex,
    std::make_index_sequence<SymmetryClass::size>{}
  );
}

//! Helper to perform rotationPeriodicity in constexpr fashion
template<typename SymmetryClass, unsigned rotationFunctionIndex>
constexpr unsigned rotationPeriodicityImpl(
  const ArrayType<unsigned, SymmetryClass::size>& runningIndices,
  const unsigned count
) {
  if(
    temple::arraysEqual(
      runningIndices,
      startingIndexSequence<SymmetryClass>()
    )
  ) {
    return count;
  }

  return rotationPeriodicityImpl<SymmetryClass, rotationFunctionIndex>(
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
constexpr unsigned rotationPeriodicity() {
  return rotationPeriodicityImpl<SymmetryClass, rotationFunctionIndex>(
    applyRotation<SymmetryClass>(
      startingIndexSequence<SymmetryClass>(),
      rotationFunctionIndex
    ),
    1
  );
}

//! Helper function to calculate all rotation periodicities
template<typename SymmetryClass, size_t ...Inds>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()>
rotationPeriodicitiesImpl(std::index_sequence<Inds...> /* inds */) {
  return {
    rotationPeriodicity<SymmetryClass, Inds>()...
  };
}

//! Calculates all multiplicities of a symmetry group's rotations.
template<typename SymmetryClass>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()> rotationPeriodicities() {
  return rotationPeriodicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
}

/*!
 * Template function generating all rotation multiplicities for a specified
 * symmetry group type
 */
template<typename SymmetryClass>
struct allRotationPeriodicities {
  static constexpr ArrayType<
    unsigned,
    SymmetryClass::rotations
  > value = rotationPeriodicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
};

//! Finds the largest size value of a set of symmetries
template<typename ... SymmetryClasses>
struct maxSymmetrySizeFunctor {
  static constexpr unsigned value() {
    ArrayType<unsigned, sizeof...(SymmetryClasses)> sizes {
      SymmetryClasses::size...
    };

    return temple::max(sizes);
  }
};

//! The largest symmetry size defined in the library
constexpr unsigned maxSymmetrySize = temple::TupleType::unpackToFunction<
  Symmetry::data::allSymmetryDataTypes,
  maxSymmetrySizeFunctor
>();

/*!
 * Fetches the coordinates of an index in a Symmetry, properly handling the
 * boost::none -> origin mapping.
 */
template<typename SymmetryClass>
constexpr temple::Vector getCoordinates(const unsigned indexInSymmetry) {
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
  const temple::Vector& i,
  const temple::Vector& j,
  const temple::Vector& k,
  const temple::Vector& l
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
    temple::Math::max(
      SymmetryClassFrom::size,
      SymmetryClassTo::size
    )
  >& indexMapping
) {
  double angularDistortion = 0;

  for(unsigned i = 0; i < SymmetryClassFrom::size; ++i) {
    for(unsigned j = i + 1; j < SymmetryClassFrom::size; ++j) {
      angularDistortion += temple::Math::abs(
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

/*! Calculates distortion caused by a symmetry transition mapping
 *
 * @param indexMapping An integer sequence specifying how indices from the
 *   source symmetry are mapped to the target symmetry
 * @param sourceSymmetrySize The size of the source symmetry
 * @param sourceAngleFunction A pointer to the source symmetry's angle function
 * @param targetAngleFunction A pointer to the target symmetry's angle function
 */
template<size_t size>
constexpr double calculateAngularDistortion(
  const ArrayType<unsigned, size>& indexMapping,
  const size_t sourceSymmetrySize,
  const data::AngleFunctionPtr& sourceAngleFunction,
  const data::AngleFunctionPtr& targetAngleFunction
) {
  double distortionSum = 0;

  for(unsigned i = 0; i < sourceSymmetrySize; ++i) {
    for(unsigned j = i + 1; j < sourceSymmetrySize; ++j) {
      distortionSum += temple::Math::abs(
        sourceAngleFunction(i, j)
        - targetAngleFunction(
          indexMapping.at(i),
          indexMapping.at(j)
        )
      );
    }
  }

  return distortionSum;
}

/*!
 * Propagates a source symmetry position through an index mapping, properly
 * handling the origin placeholder special unsigned value.
 *
 * @param symmetryPosition The symmetry position to be mapped
 * @param indexMapping The index mapping that specifies how indices are mapped
 *   from a source symmetry to a target symmetry
 */
template<size_t size>
constexpr unsigned propagateSymmetryPosition(
  const unsigned symmetryPosition,
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
 *
 * @param indexMapping The index mapping that specifies how indices are mapped
 *   from a source symmetry to a target symmetry
 */
template<typename SymmetryClassFrom, typename SymmetryClassTo>
constexpr double calculateChiralDistortion(
  const ArrayType<
    unsigned,
    temple::Math::max(
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

    chiralDistortion += temple::Math::abs(
      getTetrahedronVolume(
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(0)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(1)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(2)),
        getCoordinates<SymmetryClassFrom>(tetrahedron.at(3))
      ) - getTetrahedronVolume(
        getCoordinates<SymmetryClassTo>(
          propagateSymmetryPosition(tetrahedron.at(0), indexMapping)
        ),
        getCoordinates<SymmetryClassTo>(
          propagateSymmetryPosition(tetrahedron.at(1), indexMapping)
        ),
        getCoordinates<SymmetryClassTo>(
          propagateSymmetryPosition(tetrahedron.at(2), indexMapping)
        ),
        getCoordinates<SymmetryClassTo>(
          propagateSymmetryPosition(tetrahedron.at(3), indexMapping)
        )
      )
    );
  }

  return chiralDistortion;
}

/*!
 * Writes the indices of the original symmetry in the mapping into the target
 * symmetry's indexing scheme.
 *
 * @param mapping An index mapping that specifies how indices are mapped
 *   from a source symmetry to a target symmetry
 */
template<size_t size>
constexpr ArrayType<unsigned, size> symPosMapping(
  const ArrayType<unsigned, size>& mapping
) {
  /* Creates the list of indices in the target symmetry. Why is this necessary?
   *
   * E.g. An index mapping from linear to T-shaped. The individual
   * symmetry-internal numbering schemes are shown for the symmetry positions.
   *
   *  1  –▶  0
   *  |      |
   * (_)    (_) – 1 (new)
   *  |      |                Linear pos. 0 to
   *  0  –▶  2                pos. 2 in Tshaped
   *                                 |  ┌– Linear pos. 1 to pos. 0 in Tshaped
   *                                 |  |  ┌– This position is new
   * This mapping is represented as {2, 0, 1}.
   *
   * This function writes the indices of original mapping into the target
   * symmetry's indexing scheme.
   *
   * For this example, this returns {1, 2, 0}:
   *
   *  1 (at pos 0 in internal indexing scheme)
   *  |
   * (_) – 2 (etc.)
   *  |
   *  0
   *
   * The closely related mapping {0, 2, 1} yields target indices {0, 2, 1}.
   *
   * Which of these properties are related by target symmetry rotations?
   *
   *
   *     mapping       target indices
   * -----------------------------------
   *    {2, 0, 1}   =>   {1, 2, 0}
   *        ▲                ▲
   *        |                |
   *        X                | (C2 rotation)
   *        |                |
   *        ▼                ▼
   *    {0, 2, 1}   =>   {0, 2, 1}
   *
   */
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
 *
 * @tparam SymmetryClass A Symmetry class as defined in Primitives.h
 */
template<typename SymmetryClass>
constexpr unsigned maxRotations() {
  /* For each rotation in the symmetry class, figure out the multiplicity, i.e.
   * how often a rotation has to be applied to return to identity
   */
  constexpr auto symmetryRotationPeriodicities = rotationPeriodicities<SymmetryClass>();

  return temple::reduce(
    symmetryRotationPeriodicities,
    1u,
    std::multiplies<unsigned>()
  );
}

/*!
 * This is a perfect hash for small arrays containing numbers from 0-9. It
 * essentially interprets the array as a positive integer, e.g.:
 *
 *   {4, 8, 9, 0, 4} -> 48904
 *
 * Why? Because numeric comparisons in sets are probably faster than
 * lexicographical comparisons between vectors.
 */
template<typename UnsignedType, size_t size>
constexpr auto hashIndexList(const temple::Array<unsigned, size>& indexList) {
  constexpr unsigned maxDigitsStoreable = temple::Math::floor(
    temple::Math::log10(
      static_cast<double>(
        std::numeric_limits<UnsignedType>::max()
      )
    )
  );

  static_assert(
    size <= maxDigitsStoreable,
    "The unsigned type with which hashIndexList was instantiated cannot hold as"
    " many digits as the array to hash holds:"
  );

  UnsignedType hash = 0;
  unsigned tenPowers = 1;

  for(unsigned i = 0; i < size; ++i) {
    if(indexList.at(i) > 9) {
      throw "hashIndexList: Index list contains numbers greater than 9!";
    }
    hash += tenPowers * indexList.at(i);
    tenPowers *= 10;
  }

  return hash;
}

// Some helper types for use in generateAllRotations
template<typename SymmetryClass>
using IndicesList = ArrayType<unsigned, SymmetryClass::size>;

using IndexListStorageType = unsigned;

template<typename SymmetryClass>
using RotationsSetType = temple::DynamicSet<
  IndicesList<SymmetryClass>,
  maxRotations<SymmetryClass>()
>;

template<typename SymmetryClass>
using ChainStructuresArrayType = temple::DynamicArray<
  IndicesList<SymmetryClass>,
  maxRotations<SymmetryClass>() * 2 // TODO factor is entirely arbitrary
>;

template<typename SymmetryClass>
using ChainArrayType = temple::DynamicArray<
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
struct MappingsReturnType {
  static constexpr size_t maxMappingsSize = 50;

  using MappingsList = temple::DynamicSet<
    temple::DynamicArray<unsigned, maxSymmetrySize>,
    maxMappingsSize
  >;

  MappingsList mappings;
  double angularDistortion, chiralDistortion;

  constexpr MappingsReturnType()
    : angularDistortion(std::numeric_limits<double>::max()),
      chiralDistortion(std::numeric_limits<double>::max())
  {}

  constexpr MappingsReturnType(
    MappingsList&& passMappings,
    double&& passAngularDistortion,
    double&& passChiralDistortion
  ) : mappings(passMappings),
      angularDistortion(passAngularDistortion),
      chiralDistortion(passChiralDistortion)
  {}

  constexpr bool operator == (const MappingsReturnType& other) const {
    return (
      mappings == other.mappings
      && angularDistortion == other.angularDistortion
      && chiralDistortion == other.chiralDistortion
    );
  }

  constexpr bool operator < (const MappingsReturnType& other) const {
    return (
      std::tie(mappings, angularDistortion, chiralDistortion)
      < std::tie(other.mappings, other.angularDistortion, other.chiralDistortion)
    );
  }
};

/*!
 * Calculates the ideal index mappings when adding a ligand to a particular
 * symmetry or between symmetries of equal size.
 */
template<class SymmetryClassFrom, class SymmetryClassTo>
constexpr auto symmetryTransitionMappings() {
  static_assert(
    (
      SymmetryClassTo::size == SymmetryClassFrom::size + 1
      || SymmetryClassTo::size == SymmetryClassFrom::size
    ),
    "This function can handle only cases of equal or increasing symmetry size"
  );

  typename MappingsReturnType::MappingsList bestMappings;

  double lowestAngleDistortion = 100;
  double lowestChiralDistortion = 100;

  auto indexMapping = startingIndexSequence<SymmetryClassTo>();

  temple::Bitset<temple::Math::factorial(SymmetryClassTo::size)> encounteredMappings;

  temple::floating::ExpandedRelativeEqualityComparator<double> comparator {
    floatingPointEqualityTolerance
  };

  do {
    auto mapped = symPosMapping(indexMapping);
    auto storageVersion = temple::permutationIndex(mapped);

    if(!encounteredMappings.test(storageVersion)) {
      double angularDistortion = calculateAngularDistortion(
        indexMapping,
        SymmetryClassFrom::size,
        SymmetryClassFrom::angleFunction,
        SymmetryClassTo::angleFunction
      );

      double chiralDistortion = calculateChiralDistortion<
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
        comparator.isLessThan(angularDistortion, lowestAngleDistortion)
        || (
          comparator.isEqual(angularDistortion, lowestAngleDistortion)
          && comparator.isLessOrEqual(chiralDistortion, lowestChiralDistortion)
        )
      );

      bool clearExisting = (
        addMapping
        && !(
          comparator.isEqual(angularDistortion, lowestAngleDistortion)
          && comparator.isEqual(chiralDistortion, lowestChiralDistortion)
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
      auto allRotations = generateAllRotations<SymmetryClassTo>(mapped);

      for(const auto& rotation : allRotations) {
        encounteredMappings.set(
          temple::permutationIndex(rotation)
        );
      }
    }
  } while(temple::inPlaceNextPermutation(indexMapping));

  return MappingsReturnType(
    std::move(bestMappings),
    std::move(lowestAngleDistortion),
    std::move(lowestChiralDistortion)
  );
}

//! Find index mappings for ligand loss situations
template<class SymmetryClassFrom, class SymmetryClassTo>
constexpr auto ligandLossMappings(const unsigned deletedSymmetryPosition) {

  static_assert(
    SymmetryClassFrom::size == SymmetryClassTo::size + 1,
    "Ligand loss pathway calculation must involve symmetry size decrease"
  );

  assert(deletedSymmetryPosition < SymmetryClassFrom::size);


  typename MappingsReturnType::MappingsList bestMappings;

  double lowestAngleDistortion = 100;
  double lowestChiralDistortion = 100;

  // Construct the initial index mapping
  ArrayType<unsigned, SymmetryClassTo::size> indexMapping;
  for(unsigned i = 0; i < deletedSymmetryPosition; ++i) {
    indexMapping.at(i) = i;
  }
  for(unsigned i = deletedSymmetryPosition; i < SymmetryClassFrom::size - 1; ++i) {
    indexMapping.at(i) = i + 1;
  }

  /* NOTE: From here the algorithm is identical to symmetryTransitionMappings
   * save that symmetryTo and symmetryFrom are swapped in all occasions
   * and that inPlaceNextPermutation is only called on the subset excluding the
   * last position (the one that is added / deleted).
   */

  temple::DynamicSet<
    IndexListStorageType,
    temple::Math::factorial(SymmetryClassFrom::size)
  > encounteredMappings;

  temple::floating::ExpandedRelativeEqualityComparator<double> comparator {
    floatingPointEqualityTolerance
  };

  do {
    if(!encounteredMappings.contains(hashIndexList<IndexListStorageType>(indexMapping))) {
      double angularDistortion = calculateAngleDistortion<
        SymmetryClassTo,
        SymmetryClassFrom
      >(indexMapping);

      double chiralDistortion = calculateChiralDistortion<
        SymmetryClassTo,
        SymmetryClassFrom
      >(indexMapping);

      bool addMapping = (
        comparator.isLessThan(angularDistortion, lowestAngleDistortion)
        || (
          comparator.isEqual(angularDistortion, lowestAngleDistortion)
          && comparator.isLessOrEqual(chiralDistortion, lowestChiralDistortion)
        )
      );

      bool clearExisting = (
        addMapping
        && !(
          comparator.isEqual(angularDistortion, lowestAngleDistortion)
          && comparator.isEqual(chiralDistortion, lowestChiralDistortion)
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
      auto allRotations = generateAllRotations<SymmetryClassTo>(indexMapping);

      for(const auto& rotation : allRotations) {
        auto hashed = hashIndexList<IndexListStorageType>(rotation);

        if(!encounteredMappings.contains(hashed)) {
          encounteredMappings.insert(hashed);
        }
      }
    }
  } while(
    temple::inPlaceNextPermutation(indexMapping)
  );

  return MappingsReturnType(
    std::move(bestMappings),
    std::move(lowestAngleDistortion),
    std::move(lowestChiralDistortion)
  );
}

/* Pre-compute all ligand gain situations */
//! If symmetries are adjacent, calculate their symmetry transition mapping
template<typename SymmetrySource, typename SymmetryTarget>
constexpr
std::enable_if_t<
  (
    SymmetrySource::size == SymmetryTarget::size
    || SymmetrySource::size + 1 == SymmetryTarget::size
  ),
  temple::Optional<MappingsReturnType>
> calculateMapping() {
  return temple::Optional<MappingsReturnType> {
    symmetryTransitionMappings<SymmetrySource, SymmetryTarget>()
  };
}

//! If symmetries are not adjacent, return a None
template<typename SymmetrySource, typename SymmetryTarget>
constexpr
std::enable_if_t<
  !(
    SymmetrySource::size == SymmetryTarget::size
    || SymmetrySource::size + 1 == SymmetryTarget::size
  ),
  temple::Optional<MappingsReturnType>
> calculateMapping() {
  return temple::Optional<MappingsReturnType> {};
}

/*! Calculate stereopermutations for an unlinked symmetry
 *
 * @tparam Symmetry The symmetry for which to calculate this property.
 * @param nIdenticalLigands The number of ligands whose ranking is identical.
 *   E.g. 0 generates ABCDEF, 3 generates AAABCD, etc. for octahedral.
 */
template<typename Symmetry>
constexpr unsigned numUnlinkedAssignments(
  const unsigned nIdenticalLigands
) {
  unsigned count = 1;

  auto indices = startingIndexSequence<Symmetry>();

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  temple::DynamicSet<
    unsigned,
    temple::Math::factorial(Symmetry::size)
  > rotations;

  auto initialRotations = generateAllRotations<Symmetry>(indices);

  for(const auto& rotation : initialRotations) {
    rotations.insert(hashIndexList<unsigned>(rotation));
  }

  while(temple::inPlaceNextPermutation(indices)) {
    if(!rotations.contains(hashIndexList<unsigned>(indices))) {
      auto allRotations = generateAllRotations<Symmetry>(indices);
      for(const auto& rotation : allRotations) {
        rotations.insert(hashIndexList<unsigned>(rotation));
      }

      ++count;
    }
  }

  return count;
}

/*! Calculates whether a symmetry has multiple stereopermutations
 *
 * @tparam Symmetry The symmetry for which to calculate this property.
 * @param nIdenticalLigands The number of ligands whose ranking is identical.
 *   E.g. 0 generates ABCDEF, 3 generates AAABCD, etc. for octahedral.
 */
template<typename Symmetry>
constexpr bool hasMultipleUnlinkedAssignments(
  const unsigned nIdenticalLigands
) {
  auto indices = startingIndexSequence<Symmetry>();

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  temple::DynamicSet<
    unsigned,
    temple::Math::min(
      100u,
      temple::Math::factorial(Symmetry::size)
    )
  > rotations;

  auto initialRotations = generateAllRotations<Symmetry>(indices);

  for(const auto& rotation : initialRotations) {
    rotations.insert(hashIndexList<unsigned>(rotation));
  }

  while(temple::inPlaceNextPermutation(indices)) {
    if(!rotations.contains(hashIndexList<unsigned>(indices))) {
      return true;
    }
  }

  return false;
}

} // namespace constexprProperties

} // namespace Symmetry

#endif
