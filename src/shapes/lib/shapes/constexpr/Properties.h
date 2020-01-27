/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Compile-time shape property calculations
 *
 * Constexpr parallel to Properties.h. Contains a slew of computations on the
 * base shape data to compute derived properties. You can e.g. apply
 * rotations, extract rotation multiplicities, calculate angular and chiral
 * distortions between pairs of shapes, etc.
 *
 * @todo I don't think the constexpr ligand loss mappings are instantiated or
 * tested ANYWHERE
 */

#ifndef INCLUDE_SHAPE_CONSTEXPR_PROPERTIES_H
#define INCLUDE_SHAPE_CONSTEXPR_PROPERTIES_H

#include "shapes/Data.h"

#include "temple/constexpr/Bitset.h"
#include "temple/constexpr/DynamicSet.h"
#include "temple/constexpr/FloatingPointComparison.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Cache.h"

namespace Scine {
namespace shapes {

/**
 * @brief Compile-time calculation of shape classes' properties
 *
 * Here we calculate various properties of shape classes (i.e. classes that
 * fulfill the concepts::ShapeClass concept). Also, transition mappings can
 * be calculated at compile time using functions in this namespace.
 *
 * Data from Primitives.h takes two paths through this library: Shape classes,
 * whose equivalent data members may have different type signatures, can be
 * processed directly here and transformed. The data available in the shape
 * classes is made available to runtime functions by smoothing over the type
 * signatures, which are merely differences in array lengths, into static
 * runtime datatypes such as vectors. Then, runtime functions can calculate
 * the same properties as these available here, but by referring to shapes
 * by their Name, not by the shape class type itself.
 */
namespace constexpr_properties {

constexpr double floatingPointEqualityTolerance = 1e-4;

/*! @brief Calculates the minimum angle returned in a shape class
 *
 * @tparam ShapeClass A shape class
 *
 * @complexity{@math{\Theta(S^2)}}
 */
template<typename ShapeClass>
constexpr double calculateSmallestAngle() {
  double smallestAngle = M_PI;

  for(unsigned i = 0; i < ShapeClass::size; ++i) {
    for(unsigned j = i + 1; j < ShapeClass::size; ++j) {
      double returnedAngle = ShapeClass::angleFunction(i, j);
      if(returnedAngle < smallestAngle) {
        smallestAngle = returnedAngle;
      }
    }
  }

  return smallestAngle;
}

/**
 * @brief Metafunction calculating the smallest and largest angle that exist in
 *   a shape
 *
 * @tparam ShapeClass A class fulfilling concepts::ShapeClass
 */
template<typename ShapeClass>
struct AngleBoundsFunctor {
  /*! @brief Smallest and largest angles of the @p ShapeClass
   *
   * @complexity{@math{\Theta(S^2)}}
   */
  static constexpr std::pair<double, double> value() {
    double smallestAngle = M_PI;
    double largestAngle = 0;

    for(unsigned i = 0; i < ShapeClass::size; ++i) {
      for(unsigned j = i + 1; j < ShapeClass::size; ++j) {
        double returnedAngle = ShapeClass::angleFunction(i, j);
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
 * @brief Functor to find out the minimum angle among all shape class
 *   types passed as template arguments
 *
 * @tparam ShapeClass pack of classes fulfilling concepts::ShapeClass
 */
template<typename ... ShapeClass>
struct minAngleFunctor {
  /*! @brief Minimum angle among all shape classes
   *
   * @complexity{@math{\Theta(N S^2)}} where @math{N} is the number of shape classes and @math{S} the largest shape size
   */
  static constexpr double value() {
    const std::array<double, sizeof...(ShapeClass)> smallestAngles {{
      calculateSmallestAngle<ShapeClass>()...
    }};

    // C++17 min_element (isn't constexpr before)
    double minElement = smallestAngles.at(0);

    for(unsigned i = 1; i < sizeof...(ShapeClass); ++i) {
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
template<typename ShapeClass>
constexpr auto startingIndexSequence() {
  return temple::iota<ArrayType, unsigned, ShapeClass::size>();
}

//! Helper to perform applyRotation in constexpr fashion
template<typename ShapeClass, size_t... Indices>
constexpr ArrayType<unsigned, ShapeClass::size> applyRotationImpl(
  const ArrayType<unsigned, ShapeClass::size>& indices,
  const unsigned rotationFunctionIndex,
  std::index_sequence<Indices...> /* inds */
) {
  return {
    indices.at(
      ShapeClass::rotations.at(rotationFunctionIndex).at(Indices)
    )...
  };
}

/*! @brief Applies a shape group rotation to an array of indices
 *
 * @complexity{@math{\Theta(S)}}
 */
template<typename ShapeClass>
constexpr ArrayType<unsigned, ShapeClass::size> applyRotation(
  const ArrayType<unsigned, ShapeClass::size>& indices,
  const unsigned rotationFunctionIndex
) {
  return applyRotationImpl<ShapeClass>(
    indices,
    rotationFunctionIndex,
    std::make_index_sequence<ShapeClass::size>{}
  );
}

//! Helper to perform rotationPeriodicity in constexpr fashion
template<typename ShapeClass, unsigned rotationFunctionIndex>
constexpr unsigned rotationPeriodicityImpl(
  const ArrayType<unsigned, ShapeClass::size>& runningIndices,
  const unsigned count
) {
  if(
    temple::arraysEqual(
      runningIndices,
      startingIndexSequence<ShapeClass>()
    )
  ) {
    return count;
  }

  return rotationPeriodicityImpl<ShapeClass, rotationFunctionIndex>(
    applyRotation<ShapeClass>(
      runningIndices,
      rotationFunctionIndex
    ),
    count + 1
  );
}

/*! @brief Determines the multiplicity of a shape group rotation
 *
 * Calculates the multiplicity of a shape group's rotation specified via an
 * index in that shape's list of rotations
 *
 * @complexity{@math{\Theta(M S)} where @math{M} is the multiplicity of the
 * rotation and @math{S} is the shape size}
 *
 * @tparam ShapeClass a model of concepts::ShapeClass
 */
template<typename ShapeClass, unsigned rotationFunctionIndex>
constexpr unsigned rotationPeriodicity() {
  return rotationPeriodicityImpl<ShapeClass, rotationFunctionIndex>(
    applyRotation<ShapeClass>(
      startingIndexSequence<ShapeClass>(),
      rotationFunctionIndex
    ),
    1
  );
}

//! Helper function to calculate all rotation periodicities
template<typename ShapeClass, size_t ...Inds>
constexpr ArrayType<unsigned, ShapeClass::rotations.size()>
rotationPeriodicitiesImpl(std::index_sequence<Inds...> /* inds */) {
  return { rotationPeriodicity<ShapeClass, Inds>()... };
}

//! Calculates all multiplicities of a shape group's rotations.
template<typename ShapeClass>
constexpr ArrayType<unsigned, ShapeClass::rotations.size()> rotationPeriodicities() {
  return rotationPeriodicitiesImpl<ShapeClass>(
    std::make_index_sequence<ShapeClass::rotations.size()>{}
  );
}

/*! @brief Calculate the multiplicities of all rotations of a shape
 *
 * Template metafunction generating all rotation multiplicities for a specified
 * shape group type
 *
 * @complexity{@math{\Theta(R M S)} where @math{R} is the number of rotations
 * of the shape, @math{M} is the largest multiplicity of those rotations and
 * @math{S} is the size of the shape}
 *
 * @tparam ShapeClass a model of concepts::ShapeClass
 */
template<typename ShapeClass>
struct allRotationPeriodicities {
  static constexpr ArrayType<
    unsigned,
    ShapeClass::rotations
  > value = rotationPeriodicitiesImpl<ShapeClass>(
    std::make_index_sequence<ShapeClass::rotations.size()>{}
  );
};

//! Finds the largest size value of a set of symmetries
template<typename ... ShapeClasses>
struct maxShapeSizeFunctor {
  static constexpr unsigned value() {
    ArrayType<unsigned, sizeof...(ShapeClasses)> sizes {
      ShapeClasses::size...
    };

    return temple::max(sizes);
  }
};

//! The largest shape size defined in the library
constexpr unsigned maxShapeSize = temple::tuples::unpackToFunction<
  shapes::data::allShapeDataTypes,
  maxShapeSizeFunctor
>();

/*! @brief Fetches the coordinates of an index in a shape
 *
 * Fetches the coordinates of an index in a shape, properly handling the
 * boost::none -> origin mapping.
 *
 * @complexity{@math{\Theta(1)}}
 */
template<typename ShapeClass>
constexpr temple::Vector getCoordinates(const unsigned indexInSymmetry) {
  if(indexInSymmetry != ORIGIN_PLACEHOLDER) {
    return ShapeClass::coordinates.at(indexInSymmetry);
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

// Typedef to avoid reusing C-Style function ptr type
using AngleFunctionPtr = double(*)(const unsigned, const unsigned);

/*! @brief Calculates angular distortion caused by a shape transition mapping
 *
 * @complexity{@math{\Theta(S^2)}}
 *
 * @param indexMapping An integer sequence specifying how indices from the
 *   source shape are mapped to the target shape
 * @param sourceSymmetrySize The size of the source shape
 * @param sourceAngleFunction A pointer to the source shape's angle function
 * @param targetAngleFunction A pointer to the target shape's angle function
 */
template<size_t size>
constexpr double calculateAngularDistortion(
  const ArrayType<unsigned, size>& indexMapping,
  const size_t sourceSymmetrySize,
  const AngleFunctionPtr sourceAngleFunction,
  const AngleFunctionPtr targetAngleFunction
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

/*! @brief Propagates shape positions trhough an index mapping
 *
 * Propagates a source shape position through an index mapping, properly
 * handling the origin placeholder special unsigned value.
 *
 * @complexity{@math{\Theta(S)}}
 *
 * @param symmetryPosition The shape position to be mapped
 * @param indexMapping The index mapping that specifies how indices are mapped
 *   from a source shape to a target shape
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

/*! @brief Calculates the chiral distortion caused by a shape transition
 *
 * Calculate the chiral distortion between the source shape and the target
 * shape specified by an index mapping.
 *
 * @complexity{@math{\Theta(T)} where @math{T} is the number of tetrahedra for
 * the shape, typically small}
 *
 * @param indexMapping The index mapping that specifies how indices are mapped
 *   from a source shape to a target shape
 */
template<typename ShapeClassFrom, typename ShapeClassTo>
constexpr double calculateChiralDistortion(
  const ArrayType<
    unsigned,
    temple::Math::max(
      ShapeClassFrom::size,
      ShapeClassTo::size
    )
  >& indexMapping
) {
  double chiralDistortion = 0;

  // C++17:
  // for(const auto& tetrahedron : ShapeClassFrom::tetrahedra) {
  for(unsigned i = 0; i < ShapeClassFrom::tetrahedra.size(); ++i) {
    const auto& tetrahedron = ShapeClassFrom::tetrahedra.at(i);

    chiralDistortion += temple::Math::abs(
      getTetrahedronVolume(
        getCoordinates<ShapeClassFrom>(tetrahedron.at(0)),
        getCoordinates<ShapeClassFrom>(tetrahedron.at(1)),
        getCoordinates<ShapeClassFrom>(tetrahedron.at(2)),
        getCoordinates<ShapeClassFrom>(tetrahedron.at(3))
      ) - getTetrahedronVolume(
        getCoordinates<ShapeClassTo>(
          propagateSymmetryPosition(tetrahedron.at(0), indexMapping)
        ),
        getCoordinates<ShapeClassTo>(
          propagateSymmetryPosition(tetrahedron.at(1), indexMapping)
        ),
        getCoordinates<ShapeClassTo>(
          propagateSymmetryPosition(tetrahedron.at(2), indexMapping)
        ),
        getCoordinates<ShapeClassTo>(
          propagateSymmetryPosition(tetrahedron.at(3), indexMapping)
        )
      )
    );
  }

  return chiralDistortion;
}

/*! @brief Transform shape positions through a mapping
 *
 * Writes the indices of the original shape in the mapping into the target
 * shape's indexing scheme.
 *
 * @complexity{@math{\Theta(S)}}
 *
 * @param mapping An index mapping that specifies how indices are mapped
 *   from a source shape to a target shape
 */
template<size_t size>
constexpr ArrayType<unsigned, size> symPosMapping(
  const ArrayType<unsigned, size>& mapping
) {
  /* Creates the list of indices in the target shape. Why is this necessary?
   *
   * E.g. An index mapping from linear to T-shaped. The individual
   * shape-internal numbering schemes are shown for the shape positions.
   *
   *  1  –▶  0
   *  |      |
   * (_)    (_) – 1 (new)
   *  |      |                Line pos. 0 to
   *  0  –▶  2                pos. 2 in Tshaped
   *                                 |  ┌– Line pos. 1 to pos. 0 in Tshaped
   *                                 |  |  ┌– This position is new
   * This mapping is represented as {2, 0, 1}.
   *
   * This function writes the indices of original mapping into the target
   * shape's indexing scheme.
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
   * Which of these properties are related by target shape rotations?
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
 * Calculate a lower bound on non-superimposable rotations a shape
 * class can produce for entirely unequal indices
 *
 * @complexity{@math{\Theta(R M S)} where @math{R} is the number of rotations
 * of the shape, @math{M} is the largest multiplicity of those rotations and
 * @math{S} is the size of the shape}
 *
 * @tparam ShapeClass A shape class as defined in Primitives.h
 */
template<typename ShapeClass>
constexpr unsigned maxRotations() {
  return temple::reduce(rotationPeriodicities<ShapeClass>(), 1u, std::multiplies<>());
}

// Some helper types for use in generateAllRotations
template<typename ShapeClass>
using IndicesList = ArrayType<unsigned, ShapeClass::size>;

using IndexListStorageType = unsigned;

template<typename ShapeClass>
using RotationsSetType = temple::DynamicSet<
  IndicesList<ShapeClass>,
  maxRotations<ShapeClass>() * 2
>;

template<typename ShapeClass>
using ChainStructuresArrayType = temple::DynamicArray<
  IndicesList<ShapeClass>,
  maxRotations<ShapeClass>() * 2 // factor is entirely arbitrary
>;

template<typename ShapeClass>
using ChainArrayType = temple::DynamicArray<
  unsigned,
  maxRotations<ShapeClass>() * 2 // factor is entirely arbitrary
>;

/*! @brief Generates all rotations of a sequence of indices within a shape group
 *
 * @tparam ShapeClass a model of concepts::ShapeClass
 *
 * @complexity{At most maxRotation iterations}
 */
template<typename ShapeClass>
constexpr auto generateAllRotations(const IndicesList<ShapeClass>& indices) {
  RotationsSetType<ShapeClass> rotations;

  rotations.insert(indices);

  ChainStructuresArrayType<ShapeClass> chainStructures;
  chainStructures.push_back(indices);

  ChainArrayType<ShapeClass> chain {0u};

  while(chain.front() < ShapeClass::rotations.size()) {

    auto generated = applyRotation<ShapeClass>(
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
        && chain.back() == ShapeClass::rotations.size() - 1
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
    temple::DynamicArray<unsigned, maxShapeSize>,
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
 * added or the shape size stays the same.
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @tparam ShapeClassFrom A model of concepts::ShapeClass
 * @tparam ShapeClassTo A model of concepts::ShapeClass
 */
template<class ShapeClassFrom, class ShapeClassTo>
constexpr auto symmetryTransitionMappings() {
  static_assert(
    (
      ShapeClassTo::size == ShapeClassFrom::size + 1
      || ShapeClassTo::size == ShapeClassFrom::size
    ),
    "This function can handle only cases of equal or increasing shape size"
  );

  typename MappingsReturnType::MappingsList bestMappings;

  double lowestAngleDistortion = 100;
  double lowestChiralDistortion = 100;

  auto indexMapping = startingIndexSequence<ShapeClassTo>();

  temple::Bitset<temple::Math::factorial(ShapeClassTo::size)> encounteredMappings;

  temple::floating::ExpandedRelativeEqualityComparator<double> comparator {
    floatingPointEqualityTolerance
  };

  do {
    auto mapped = symPosMapping(indexMapping);
    auto storageVersion = temple::permutationIndex(mapped);

    if(!encounteredMappings.test(storageVersion)) {
      double angularDistortion = calculateAngularDistortion(
        indexMapping,
        ShapeClassFrom::size,
        ShapeClassFrom::angleFunction,
        ShapeClassTo::angleFunction
      );

      double chiralDistortion = calculateChiralDistortion<
        ShapeClassFrom,
        ShapeClassTo
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
      auto allRotations = generateAllRotations<ShapeClassTo>(mapped);

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
 * @tparam ShapeClassFrom A model of concepts::ShapeClass
 * @tparam ShapeClassTo A model of concepts::ShapeClass
 */
template<class ShapeClassFrom, class ShapeClassTo>
constexpr auto ligandLossMappings(const unsigned deletedSymmetryPosition) {

  static_assert(
    ShapeClassFrom::size == ShapeClassTo::size + 1,
    "Ligand loss pathway calculation must involve shape size decrease"
  );

  assert(deletedSymmetryPosition < ShapeClassFrom::size);


  typename MappingsReturnType::MappingsList bestMappings;

  double lowestAngleDistortion = 100;
  double lowestChiralDistortion = 100;

  // Construct the initial index mapping
  ArrayType<unsigned, ShapeClassTo::size> indexMapping;
  for(unsigned i = 0; i < deletedSymmetryPosition; ++i) {
    indexMapping.at(i) = i;
  }
  for(unsigned i = deletedSymmetryPosition; i < ShapeClassFrom::size - 1; ++i) {
    indexMapping.at(i) = i + 1;
  }

  /* NOTE: From here the algorithm is identical to symmetryTransitionMappings
   * save that symmetryTo and symmetryFrom are swapped in all occasions
   */
  temple::Bitset<temple::Math::factorial(ShapeClassFrom::size)> encounteredMappings;

  temple::floating::ExpandedRelativeEqualityComparator<double> comparator {
    floatingPointEqualityTolerance
  };

  do {
    auto mapped = symPosMapping(indexMapping);
    auto storageVersion = temple::permutationIndex(mapped);

    if(!encounteredMappings.test(storageVersion)) {
      double angularDistortion = calculateAngularDistortion(
        indexMapping,
        ShapeClassTo::size,
        ShapeClassTo::angleFunction,
        ShapeClassFrom::angleFunction
      );

      double chiralDistortion = calculateChiralDistortion<
        ShapeClassTo,
        ShapeClassFrom
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
      auto allRotations = generateAllRotations<ShapeClassTo>(indexMapping);

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
//! If symmetries are adjacent, calculate their shape transition mapping
template<typename SymmetrySource, typename SymmetryTarget>
constexpr
std::enable_if_t<
  (
    SymmetryTarget::size <= 7
    && (
      SymmetrySource::size == SymmetryTarget::size
      || SymmetrySource::size + 1 == SymmetryTarget::size
    )
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
    SymmetryTarget::size <= 7
    && (
      SymmetrySource::size == SymmetryTarget::size
      || SymmetrySource::size + 1 == SymmetryTarget::size
    )
  ),
  temple::Optional<MappingsReturnType>
> calculateMapping() {
  return temple::Optional<MappingsReturnType> {};
}

/*!
 * @brief Calculate stereopermutations for an unlinked shape
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @tparam ShapeClass A model of concepts::ShapeClass
 * @param nIdenticalLigands The number of ligands whose ranking is identical.
 *   E.g. 0 generates ABCDEF, 3 generates AAABCD, etc. for octahedral.
 */
template<typename ShapeClass>
constexpr unsigned numUnlinkedStereopermutations(
  const unsigned nIdenticalLigands
) {
  unsigned count = 1;

  auto indices = startingIndexSequence<ShapeClass>();

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  temple::Bitset<temple::Math::factorial(ShapeClass::size)> rotations;

  auto initialRotations = generateAllRotations<ShapeClass>(indices);

  for(const auto& rotation : initialRotations) {
    rotations.set(temple::permutationIndex(rotation));
  }

  while(temple::inPlaceNextPermutation(indices)) {
    if(!rotations.test(temple::permutationIndex(indices))) {
      auto allRotations = generateAllRotations<ShapeClass>(indices);
      for(const auto& rotation : allRotations) {
        rotations.set(temple::permutationIndex(rotation));
      }

      ++count;
    }
  }

  return count;
}

/*! @brief Calculates whether a shape has multiple stereopermutations
 *
 * @complexity{@math{\Theta(S!)}}
 *
 * @tparam ShapeClass The shape for which to calculate this property.
 * @param nIdenticalLigands The number of ligands whose ranking is identical.
 *   E.g. 0 generates ABCDEF, 3 generates AAABCD, etc. for octahedral.
 */
template<typename ShapeClass>
constexpr bool hasMultipleUnlinkedStereopermutations(const unsigned nIdenticalLigands) {
  if(nIdenticalLigands == ShapeClass::size) {
    return false;
  }

  auto indices = startingIndexSequence<ShapeClass>();

  for(unsigned i = 0; i < nIdenticalLigands; ++i) {
    indices.at(i) = 0;
  }

  temple::Bitset<temple::Math::factorial(ShapeClass::size)> rotations;

  auto initialRotations = generateAllRotations<ShapeClass>(indices);

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

} // namespace constexpr_properties
} // namespace shapes
} // namespace Scine

#endif
