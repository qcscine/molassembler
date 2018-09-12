#ifndef INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_SYMMETRIES_H
#define INCLUDE_MOLASSEMBLER_CHEMICAL_SYMMETRIES_SYMMETRIES_H

#include "boost/optional.hpp"
#include "Eigen/Core"
#include "temple/constexpr/TupleType.h"

#include "chemical_symmetries/Primitives.h"

#include <vector>
#include <map>

/*! @file
 *
 * @brief Centralizes basic symmetry data in name-accessible fashion
 *
 * Central inclusion point of the library. Defines the main symmetry data and
 * all accessors. Symmetries are defined in a constexpr fashion and homogenized
 * into a single container at compile-time to allow for compile-time computation
 * without losing universal accessibility at run-time.
 */

namespace Symmetry {

/* Typedefs */
//! The type to store symmetry rotations
using RotationsList = std::vector<
  std::vector<unsigned>
>;

/*!
 * All angle functions can be called with arbitrary (valid) parameters
 * without failing. Valid here means that a != b and less than the size of
 * the symmetry requested.
 *
 * They return angles in radians.
 */
using AngleFunctionType = std::function<
  double(const unsigned, const unsigned)
>;

/*!
 * All symmetries have a guess implementation of what could work as the defined
 * tetrahedra. Have to use boost::none to signal to replace this position with
 * the central atom as it is not part of the indexing scheme used here.
 *
 * In case all higher symmetries than trigonal pyramidal are representable
 * without boost::none and that proves to work, then perhaps make an exception
 * for it and treat all others without the optional. If that cannot be done,
 * consider refactoring (changing the numbering scheme in some fashion that
 * boost::none does not have to be used.
 */
using TetrahedronList = std::vector<
  std::array<
    boost::optional<unsigned>,
    4
  >
>;

using CoordinateList = std::vector<Eigen::Vector3d>;

//! Dynamic symmetry information data struct
struct SymmetryInformation {
  const std::string stringName;
  const unsigned size;
  const RotationsList rotations;
  const TetrahedronList tetrahedra;
  const CoordinateList coordinates;

  // Direct initialization
  SymmetryInformation(
    std::string passStringName,
    unsigned passSize,
    RotationsList passRotations,
    TetrahedronList passTetrahedra,
    CoordinateList passCoordinates
  ) : stringName(std::move(passStringName)),
      size(passSize),
      rotations(std::move(passRotations)),
      tetrahedra(std::move(passTetrahedra)),
      coordinates(std::move(passCoordinates))
  {}
};

// Helper function to create all names vector
template<size_t ... Inds>
constexpr std::array<Name, nSymmetries> makeAllNames(
  std::index_sequence<Inds...> /* indexSequence */
) {
  return {{
    static_cast<Name>(Inds)...
  }};
}

//! A list of all the enum class values
constexpr std::array<Name, nSymmetries> allNames = makeAllNames(
  std::make_index_sequence<nSymmetries>()
);

namespace data {

// Typedef to avoid reusing C-Style function ptr type
using AngleFunctionPtr = double(*)(const unsigned, const unsigned);

/*!
 * Constructs an array of function pointers to all static angle functions
 * for runtime lookup
 */
template<typename ...SymmetryClasses>
struct angleFunctionFunctor {
  static constexpr std::array<
    data::AngleFunctionPtr,
    sizeof...(SymmetryClasses)
  > value() {
    return {{
      &SymmetryClasses::angleFunction...
    }};
  }
};

/*! Conversion function to make the dynamic rotations list type from the
 * constexpr data types given in a specifc symmetry class type
 */
template<size_t symmetrySize, size_t nRotations>
RotationsList makeRotations(
  const std::array<
    std::array<unsigned, symmetrySize>,
    nRotations
  >& constexprRotations
) {
  RotationsList rotations;

  for(const auto& rotation : constexprRotations) {
    rotations.emplace_back(
      rotation.begin(),
      rotation.end()
    );
  }

  return rotations;
}

/*! Conversion function to make the dynamic tetrahedron list type from the
 * constexpr data types given in a specifc symmetry class type
 */
template<size_t nTetrahedra>
TetrahedronList makeTetrahedra(
  const std::array<
    std::array<unsigned, 4>,
    nTetrahedra
  >& constexprTetrahedra
) {
  TetrahedronList tetrahedra (nTetrahedra);

  for(unsigned i = 0; i < nTetrahedra; ++i) {
    for(unsigned j = 0; j < 4; ++j) {
      const auto& constexprValue = constexprTetrahedra.at(i).at(j);

      if(constexprValue == ORIGIN_PLACEHOLDER) {
        tetrahedra.at(i).at(j) = boost::none;
      } else {
        tetrahedra.at(i).at(j) = constexprValue;
      }
    }
  }

  return tetrahedra;
}

//! Conversion helper to Eigen type from constexpr vector type
Eigen::Vector3d toEigen(const temple::Vector& cVector);

/*! Conversion function to make the dynamic coordinates list type from the
 * constexpr data types given in a specifc symmetry class type
 */
template<size_t symmetrySize>
CoordinateList makeCoordinates(
  const std::array<temple::Vector, symmetrySize>& constexprCoordinates
) {
  CoordinateList coordinates (symmetrySize);

  for(unsigned i = 0; i < symmetrySize; ++i) {
    coordinates.at(i) = toEigen(constexprCoordinates.at(i));
  }

  return coordinates;
}

/*! This constructs the SymmetryInformation instance from a specific symmetry
 * class type.
 */
template<typename SymmetryClass>
SymmetryInformation makeSymmetryInformation() {
  return {
    SymmetryClass::stringName,
    SymmetryClass::size,
    makeRotations(SymmetryClass::rotations),
    makeTetrahedra(SymmetryClass::tetrahedra),
    makeCoordinates(SymmetryClass::coordinates)
  };
}

/*! This creates a map initialization pair for a specific symmetry class type.
 * The key is the name, the mapped_type a SymmetryInformation instance
 */
template<typename SymmetryClass>
std::pair<Name, SymmetryInformation> makeMapInitPair() {
  return {
    SymmetryClass::name,
    makeSymmetryInformation<SymmetryClass>()
  };
}

/*! Creates the mapping between a symmetry class's name and it's dynamic
 * information in order to have runtime lookup based on symmetry names.
 */
template<typename ...SymmetryClasses>
struct symmetryInformationFunctor {
  static const std::map<Name, SymmetryInformation> value() {
    return {{
      makeMapInitPair<SymmetryClasses>()...
    }};
  }
};

//! An array containing pointers to all symmetry data types' angle function
constexpr auto angleFunctions = temple::TupleType::unpackToFunction<
  allSymmetryDataTypes,
  angleFunctionFunctor
>();

} // namespace data

/* Core symmetry data, this has dynamic types and is hence initialized in the
 * .cpp file from the tuple containing all symmetry data types and the
 * symmetryInformationFunctor
 */
const std::map<Name, SymmetryInformation>& symmetryData();

/* Interface */
//! Fetch the string name of a symmetry
inline const std::string& name(const Name name) {
  return symmetryData().at(name).stringName;
}

inline Name nameFromString(const std::string& nameString) {
  for(const Name symmetryName : allNames) {
    if(symmetryData().at(symmetryName).stringName == nameString) {
      return symmetryName;
    }
  }

  throw std::logic_error("No symmetry exists under that name!");
}

//! Fetch a space-free name for file naming
std::string spaceFreeName(Name name);

//! Fetch the number of symmetry positions of a symmetry
inline unsigned size(const Name name) {
  return symmetryData().at(name).size;
}

//! Fetches a symmetry's list of rotations
inline const RotationsList& rotations(const Name name) {
  return symmetryData().at(name).rotations;
}

//! Gets a symmetry's angle function
inline data::AngleFunctionPtr angleFunction(const Name name) {
  auto symmetryIndex = static_cast<unsigned>(name);
  return data::angleFunctions.at(symmetryIndex);
}

//! Returns the index of a symmetry name within allNames
unsigned nameIndex(Name name) PURITY_STRONG;

//! Fetches the list of tetrahedra defined in a symmetry
inline const TetrahedronList& tetrahedra(const Name name) {
  return symmetryData().at(name).tetrahedra;
}

} // namespace Symmetry

#endif
