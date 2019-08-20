/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Compile-time symmetry property calculations
 *
 * Constexpr parallel to Properties.h. Contains a slew of computations on the
 * base symmetry data to compute derived properties. You can e.g. apply
 * rotations, extract rotation multiplicities, calculate angular and chiral
 * distortions between pairs of symmetries, etc.
 *
 * @todo I don't think the constexpr ligand loss mappings are instantiated or
 * tested ANYWHERE
 */

#ifndef INCLUDE_SYMMETRY_INFORMATION_CONSTEXPR_PROPERTIES_H
#define INCLUDE_SYMMETRY_INFORMATION_CONSTEXPR_PROPERTIES_H

#include "chemical_symmetries/Symmetries.h"

#include "temple/constexpr/Bitset.h"
#include "temple/constexpr/DynamicSet.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Cache.h"

namespace Scine {

namespace Symmetry {

/**
 * @brief Compile-time amenable calculation of symmetry classes' properties
 *
 * Here we calculate various properties of symmetry classes (i.e. classes that
 * fulfill the concepts::SymmetryClass concept). Also, transition mappings can
 * be calculated at compile time using functions in this namespace.
 *
 * Data from Primitives.h takes two paths through this library: Symmetry classes,
 * whose equivalent data members may have different type signatures, can be
 * processed directly here and transformed. The data available in the symmetry
 * classes is made available to runtime functions by smoothing over the type
 * signatures, which are merely differences in array lengths, into static
 * runtime datatypes such as vectors. Then, runtime functions can calculate
 * the same properties as these available here, but by referring to symmetries
 * by their Name, not by the symmetry class type itself.
 */
namespace constexprProperties {

constexpr double floatingPointEqualityTolerance = 1e-4;

/*! @brief Calculates the minimum angle returned in a symmetry class
 *
 * @tparam SymmetryClass A Symmetry class
 *
 * @complexity{@math{\Theta(S^2)}}
 */
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

/**
 * @brief Metafunction calculating the smallest and largest angle that exist in
 *   a symmetry
 *
 * @tparam SymmetryClass A class fulfilling concepts::SymmetryClass
 */
template<typename SymmetryClass>
struct AngleBoundsFunctor {
  /*! @brief Smallest and largest angles of the @p SymmetryClass
   *
   * @complexity{@math{\Theta(S^2)}}
   */
  static constexpr std::pair<double, double> value() {
    double smallestAngle = M_PI;
    double largestAngle = 0;

    for(unsigned i = 0; i < SymmetryClass::size; ++i) {
      for(unsigned j = i + 1; j < SymmetryClass::size; ++j) {
        double returnedAngle = SymmetryClass::angleFunction(i, j);
        if(returnedAngle < smallestAngle) {
          smallestAngle = returnedAngle;
        }

        if(returnedAngle > largestAngle) {
          largestAngle = returnedAngle;
        }
      }
    }

    return {smallestAngle, largestAngle};
  }
};

/*!
 * @brief Functor to find out the minimum angle among all symmetry class
 *   types passed as template arguments
 *
 * @tparam SymmetryClasses pack of classes fulfilling concepts::SymmetryClass
 */
template<typename ... SymmetryClasses>
struct minAngleFunctor {
  /*! @brief Minimum angle among all symmetry classes
   *
   * @complexity{@math{\Theta(N S^2)}} where @math{N} is the number of symmetry classes and @math{S} the largest symmetry size
   */
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
  return temple::iota<ArrayType, unsigned, SymmetryClass::size>();
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

/*! @brief Applies a symmetry group rotation to an array of indices
 *
 * @complexity{@math{\Theta(S)}}
 */
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

/*! @brief Determines the multiplicity of a symmetry group rotation
 *
 * Calculates the multiplicity of a symmetry group's rotation specified via an
 * index in that symmetry's list of rotations
 *
 * @complexity{@math{\Theta(M S)} where @math{M} is the multiplicity of the
 * rotation and @math{S} is the symmetry size}
 *
 * @tparam SymmetryClass a model of concepts::SymmetryClass
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
  return { rotationPeriodicity<SymmetryClass, Inds>()... };
}

//! Calculates all multiplicities of a symmetry group's rotations.
template<typename SymmetryClass>
constexpr ArrayType<unsigned, SymmetryClass::rotations.size()> rotationPeriodicities() {
  return rotationPeriodicitiesImpl<SymmetryClass>(
    std::make_index_sequence<SymmetryClass::rotations.size()>{}
  );
}

/*! @brief Calculate the multiplicities of all rotations of a symmetry
 *
 * Template metafunction generating all rotation multiplicities for a specified
 * symmetry group type
 *
 * @complexity{@math{\Theta(R M S)} where @math{R} is the number of rotations
 * of the symmetry, @math{M} is the largest multiplicity of those rotations and
 * @math{S} is the size of the symmetry}
 *
 * @tparam SymmetryClass a model of concepts::SymmetryClass
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

/*! @brief Fetches the coordinates of an index in a Symmetry
 *
 * Fetches the coordinates of an index in a Symmetry, properly handling the
 * boost::none -> origin mapping.
 *
 * @complexity{@math{\Theta(1)}}
 */
template<typename SymmetryClass>
constexpr temple::Vector getCoordinates(const unsigned indexInSymmetry) {
  if(indexInSymmetry != ORIGIN_PLACEHOLDER) {
    return SymmetryClass::coordinates.at(indexInSymmetry);
  }

  return {0, 0, 0};
}

/*! @brief Calculates the volume of a tetrahedron spanned by four positions
 *
 * @complexity{@math{\Theta(1)}}
 *
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

/*! @brief Calculates angular distortion caused by a symmetry transition mapping
 *
 * @complexity{@math{\Theta(S^2)}}
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

/*! @brief Propagates symmetry positions trhough an index mapping
 *
 * Propagates a source symmetry position through an index mapping, properly
 * handling the origin placeholder special unsigned value.
 *
 * @complexity{@math{\Theta(S)}}
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

/*! @brief Calculates the chiral distortion caused by a symmetry transition
 *
 * Calculate the chiral distortion between the source symmetry and the target
 * symmetry specified by an index mapping.
 *
 * @complexity{@math{\Theta(T)} where @math{T} is the number of tetrahedra for
 * the symmetry, typically small}
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

/*! @brief Transform symmetry positions through a mapping
 *
 * Writes the indices of the original symmetry in the mapping into the target
 * symmetry's indexing scheme.
 *
 * @complexity{@math{\Theta(S)}}
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
 * Calculate an upper bound on non-superimposable rotations a symmetry
 * class can produce for entirely unequal indices
 *
 * @complexity{@math{\Theta(R M S)} where @math{R} is the number of rotations
 * of the symmetry, @math{M} is the largest multiplicity of those rotations and
 * @math{S} is the size of the symmetry}
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
    std::multiplies<>()
  );
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
  maxRotations<SymmetryClass>() * 2 // factor is entirely arbitrary
>;

template<typename SymmetryClass>
using ChainArrayType = temple::DynamicArray<
  unsigned,
  maxRotations<SymmetryClass>() * 2 // factor is entirely arbitrary
>;

/*! @brief Generates all rotations of a sequence of indices within a symmetry group
 *
 * @tparam SymmetryClass a model of concepts::SymmetryClass
 *
 * @complexity{At most maxRotation iterations}
 */
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
 * @brief Data struct to collect the results of calculating the ideal index
 *   mappings between pairs of indices
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

  //! Lexicographical equivalence
  constexpr bool operator == (const MappingsReturnType& other) const {
    return (
      mappings == other.mappings
      && angularDistortion == other.angularDistortion
      && chiralDistortion == other.chiralDistortion
    );
  }

  //! Lexicographical ordering
  constexpr bool operator < (const MappingsReturnType& other) const {
    return (
      std::tie(mappings, angularDistortion, chiralDistortion)
      < std::tie(other.mappings, other.angularDistortion, other.chiralDistortion)
    );
  }
};

/*! @brief Calculates ideal index mappings for +1, 0 size transitions
 *
 * Calculates the ideal index mappings for transitions in which a ligand is
 * added or the symmetry size stays the same.
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @tparam SymmetryClassFrom A model of concepts::SymmetryClass
 * @tparam SymmetryClassTo A model of concepts::SymmetryClass
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

/*! @brief Find index mappings for ligand loss situations
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @tparam SymmetryClassFrom A model of concepts::SymmetryClass
 * @tparam SymmetryClassTo A model of concepts::SymmetryClass
 */
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
   */
  temple::Bitset<temple::Math::factorial(SymmetryClassFrom::size)> encounteredMappings;

  temple::floating::ExpandedRelativeEqualityComparator<double> comparator {
    floatingPointEqualityTolerance
  };

  do {
    auto mapped = symPosMapping(indexMapping);
    auto storageVersion = temple::permutationIndex(mapped);

    if(!encounteredMappings.test(storageVersion)) {
      double angularDistortion = calculateAngularDistortion(
        indexMapping,
        SymmetryClassTo::size,
        SymmetryClassTo::angleFunction,
        SymmetryClassFrom::angleFunction
      );

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

/*!
 * @brief Calculate stereopermutations for an unlinked symmetry
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @tparam SymmetryClass A model of concepts::SymmetryClass
 * @param nIdenticalLigands The number of ligands whose ranking is identical.
 *   E.g. 0 generates ABCDEF, 3 generates AAABCD, etc. for octahedral.
 */
template<typename SymmetryClass>
constexpr unsigned numUnlinkedStereopermutations(
  const unsigned nIdenticalLigands
) {
  unsigned count = 1;

  auto indices = startingIndexSequence<SymmetryClass>();

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  temple::Bitset<temple::Math::factorial(SymmetryClass::size)> rotations;

  auto initialRotations = generateAllRotations<SymmetryClass>(indices);

  for(const auto& rotation : initialRotations) {
    rotations.set(temple::permutationIndex(rotation));
  }

  while(temple::inPlaceNextPermutation(indices)) {
    if(!rotations.test(temple::permutationIndex(indices))) {
      auto allRotations = generateAllRotations<SymmetryClass>(indices);
      for(const auto& rotation : allRotations) {
        rotations.set(temple::permutationIndex(rotation));
      }

      ++count;
    }
  }

  return count;
}

/*! @brief Calculates whether a symmetry has multiple stereopermutations
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @tparam Symmetry The symmetry for which to calculate this property.
 * @param nIdenticalLigands The number of ligands whose ranking is identical.
 *   E.g. 0 generates ABCDEF, 3 generates AAABCD, etc. for octahedral.
 */
template<typename SymmetryClass>
constexpr bool hasMultipleUnlinkedStereopermutations(const unsigned nIdenticalLigands) {
  if(nIdenticalLigands == SymmetryClass::size) {
    return false;
  }

  auto indices = startingIndexSequence<SymmetryClass>();

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  temple::Bitset<temple::Math::factorial(SymmetryClass::size)> rotations;

  auto initialRotations = generateAllRotations<SymmetryClass>(indices);

  for(const auto& rotation : initialRotations) {
    rotations.set(temple::permutationIndex(rotation));
  }

  while(temple::inPlaceNextPermutation(indices)) {
    if(!rotations.test(temple::permutationIndex(indices))) {
      return true;
    }
  }

  return false;
}

} // namespace constexprProperties

} // namespace Symmetry

} // namespace Scine

#endif
