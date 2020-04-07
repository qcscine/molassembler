/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Interpret multiple molecules in positional information
 *
 * Contains functionality permitting the interpretation of several Molecules
 * from three-dimensional structures with or without accompanying bond orders.
 */

#ifndef INCLUDE_MOLASSEMBLER_INTERPRET_H
#define INCLUDE_MOLASSEMBLER_INTERPRET_H

#include "molassembler/Export.h"
#include "boost/optional.hpp"
#include <vector>


namespace Scine {
namespace Utils {

// External forward-declarations
enum class ElementType : unsigned;
class AtomCollection;
class BondOrderCollection;
using ElementTypeCollection = std::vector<ElementType>;
} // namespace Utils

namespace molassembler {

// Forward-declarations
class Molecule;
class Graph;
class AngstromPositions;

//! @brief Given Cartesian coordinates, construct graphs or molecules
namespace interpret {

//! @brief How floating-point bond orders are discretized into bond types
enum class MASM_EXPORT BondDiscretizationOption {
  //! @brief All bond orders >= 0.5 are considered single bonds
  Binary,
  //! @brief Bond orders are rounded to the nearest integer
  RoundToNearest
};

//! Type used to represent a map from an atom collection index to an interpreted object
using ComponentMap = std::vector<unsigned>;

//! Result type of an interpret call.
struct MASM_EXPORT MoleculesResult {
  //! The list of individual molecules found in the 3D information
  std::vector<Molecule> molecules;
  //! A map from an index within the AtomCollection to its index in molecules member
  ComponentMap componentMap;
};

/*! @brief The function that actually does all the work with the library-internal wrapper
 *
 * @complexity{@math{\Theta(M)} molecule instantiations for each connected
 * component found of at least linear complexity each}
 *
 * @param elements Element type collection
 * @param angstromWrapper Positional information in Angstrom units
 * @param bondOrders Bond orders
 * @param discretization How to discretize fractional bond orders
 * @param stereopermutatorThreshold From which fractional bond
 *   order on to try the interpretation of bond stereopermutator. If set as
 *   @p boost::none, no bond stereopermutators are interpreted.
 *
 * @throws invalid_argument If the number of particles in the element
 *   collection, angstrom wrapper or bond order collection do not match.
 *
 * @returns A list of found molecules and an index mapping to each molecule
 */
MASM_EXPORT MoleculesResult molecules(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromWrapper,
  const Utils::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereopermutatorThreshold = 1.4
);

/*! @brief Interpret a molecule from positional information only. Calculates
 *   bond orders using uffBondOrders.
 *
 * @param elements Element type collection
 * @param angstromWrapper Positional information in Angstrom units
 * @param discretization How to discretize fractional bond orders
 * @param stereopermutatorThreshold From which fractional bond
 *   order on to try the interpretation of bond stereopermutator. If set as
 *   @p boost::none, no bond stereopermutators are interpreted.
 *
 * @throws invalid_argument If the number of particles in the element
 *   collection and angstrom wrapper do not match.
 *
 * @warning Using UFF bond order calculation is often not even wrong, i.e. so
 *   bad as to be completely unusable. Prefer interpreting using supplied bond
 *   orders from a more advanced calculation.
 *
 * @returns A list of found molecules and an index mapping to each molecule
 */
MASM_EXPORT MoleculesResult molecules(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromWrapper,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereopermutatorThreshold = 1.4
);

/*!
 * @brief Interpret molecules from element types, positional information and a bond order collection.
 *
 * Attempts to interpret (possibly multiple) Molecules from element types,
 * positional information and a bond order collection. Bond orders are
 * discretized into bond types. Connected components within the space are
 * identified and individually instantiated into Molecules. The instantiation
 * behavior of BondStereopermutators in the Molecules can be limited to edges whose
 * bond order exceeds a particular value.
 *
 * @param atomCollection Element types and positional information in Bohr units.
 * @param bondOrders Fractional bond orders
 * @param discretization Decide how bond orders are discretized into bond types
 * @param stereopermutatorThreshold If specified, limits the
 *   instantiation of BondStereopermutators onto edges whose fractional bond orders
 *   exceed the provided threshold. If this is not desired, specify boost::none.
 *
 * @throws invalid_argument If the number of particles in the atom
 *   collection and bond order collection do not match.
 *
 * @note Assumes that the provided atom collection's positions are in
 *   Bohr units.
 */
MASM_EXPORT MoleculesResult molecules(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereopermutatorThreshold = 1.4
);

/*!
 * @brief Interpret molecules in 3D information.
 *
 * Attempts to interpret (possibly multiple) Molecules from element types and
 * positional information. Bond orders are calculated from atom-pairwise
 * spatial distances using UFF parameters. The bond orders are then discretized
 * into bond types. Connected components within the space are identified and
 * individually instantiated into Molecules. The instantiation behavior of
 * BondStereopermutators in the Molecules can be limited to edges whose bond order
 * exceeds a particular value.
 *
 * @param atomCollection Element types and positional information in Bohr units.
 * @param discretization Decide how bond orders are discretized into bond types
 * @param stereopermutatorThreshold If specified, limits the
 *   instantiation of BondStereopermutators onto edges whose fractional bond orders
 *   exceed the provided threshold
 *
 * @note Assumes that the provided atom collection's positions are in
 * Bohr units.
 *
 * @warning UFF parameter bond order calculation is very primitive and carries
 *   a high risk of misinterpretation
 */
MASM_EXPORT MoleculesResult molecules(
  const Utils::AtomCollection& atomCollection,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereopermutatorThreshold = 1.4
);

/**
 * @brief Splits an AtomCollection just like the interpret split the positions
 *   into parts
 *
 * @param componentMap Component map received from an interpretation
 * @param atomCollection AtomCollection to split
 *
 * @return As many atom collections as molecules from interpretation
 */
MASM_EXPORT std::vector<Utils::AtomCollection> applyInterpretationMap(
  const ComponentMap& componentMap,
  const Utils::AtomCollection& atomCollection
);

/*!
 * @brief Yields mapping from indices in components to original input indices
 *
 * @param componentMap Component mapping from an interpret call
 *
 * @return List of flat maps from indices within components to the original
 *   input index
 */
MASM_EXPORT std::vector<
  std::vector<unsigned>
> invertComponentMap(const ComponentMap& componentMap);

//! Result type of a graph interpret call
struct MASM_EXPORT GraphsResult {
  //! Individual graphs found
  std::vector<Graph> graphs;
  //! Mapping from atom collection index to graph component index
  ComponentMap componentMap;
};

/*!
 * @brief The function that actually does all the work with the library-internal wrapper
 *
 * @complexity{@math{\Theta(M)} graph instantiations for each connected
 * component found of at least linear complexity each}
 *
 * @param elements Element type collection
 * @param angstromWrapper Positional information in Angstrom units
 * @param bondOrders Bond orders
 * @param discretization How to discretize fractional bond orders
 *
 * @throws invalid_argument If the number of particles in the element
 *   collection, angstrom wrapper or bond order collection do not match.
 *
 * @returns A list of found graphs and an index mapping to each graph
 */
MASM_EXPORT GraphsResult graphs(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromWrapper,
  const Utils::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary
);

/*!
 * @brief Interpret graphs from element types, positional information and a bond order collection.
 *
 * Attempts to interpret (possibly multiple) graphs from element types,
 * positional information and a bond order collection. Bond orders are
 * discretized into bond types. Connected components within the space are
 * identified and individually instantiated into Molecules.
 *
 * @param atomCollection Element types and positional information in Bohr units.
 * @param bondOrders Fractional bond orders
 * @param discretization Decide how bond orders are discretized into bond types
 *
 * @throws invalid_argument If the number of particles in the atom
 *   collection and bond order collection do not match.
 *
 * @note Assumes that the provided atom collection's positions are in
 *   Bohr units.
 */
MASM_EXPORT GraphsResult graphs(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary
);

} // namespace interpret
} // namespace molassembler
} // namespace Scine

#endif
