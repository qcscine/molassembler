/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Interpret multiple molecules in positional information
 *
 * Contains functionality permitting the interpretation of several Molecules
 * from three-dimensional structures with or without accompanying bond orders.
 */

#ifndef INCLUDE_MOLASSEMBLER_INTERPRET_H
#define INCLUDE_MOLASSEMBLER_INTERPRET_H

#include "Molassembler/Export.h"
#include "Molassembler/PeriodicBoundaries.h"
#include "boost/optional.hpp"

namespace Scine {
namespace Utils {

// External forward-declarations
enum class ElementType : unsigned;
class AtomCollection;
class BondOrderCollection;
using ElementTypeCollection = std::vector<ElementType>;
} // namespace Utils

namespace Molassembler {

// Forward-declarations
class Molecule;
class AngstromPositions;

//! @brief Given Cartesian coordinates, construct graphs or molecules
namespace Interpret {

//! @brief How floating-point bond orders are discretized into bond types
enum class MASM_EXPORT BondDiscretizationOption {
  //! @brief All bond orders >= 0.5 are considered single bonds
  Binary,
  //! @brief Bond orders are rounded to the nearest integer
  RoundToNearest
};

//! Type used to represent a map from an atom collection index to an interpreted object
struct MASM_EXPORT ComponentMap {
  struct ComponentIndexPair {
    //! Index of component
    unsigned component;
    //! Atom index
    unsigned long atomIndex;
  };

  //! Transform index from original atom collection to component and atom index
  ComponentIndexPair apply(unsigned index) const;
  //! Transform component and atom index to index from original atom collection
  unsigned invert(const ComponentIndexPair& pair) const;


  /*!
   * @brief Yields mapping from indices in components to original input indices
   * @return List of flat maps from indices within components to the original
   *   input index
   */
  std::vector<
    std::vector<unsigned>
  > invert() const;

  /**
   * @brief Splits an AtomCollection just like the interpret split the positions
   *   into parts
   *
   * @param atomCollection AtomCollection to split
   *
   * @return As many atom collections as molecules from interpretation
   */
  std::vector<Utils::AtomCollection> apply(const Utils::AtomCollection& atomCollection) const;

  MASM_NO_EXPORT std::vector<AngstromPositions> apply(const AngstromPositions& positions) const;

  MASM_NO_EXPORT std::vector<PeriodicBoundaryDuplicates> apply(
    const std::unordered_set<unsigned>& uninterestingAtoms,
    const std::unordered_map<unsigned, unsigned>& ghostAtomMap
  ) const;

  //! Size of the map
  inline unsigned size() const { return map.size(); }

  unsigned countComponents() const;

//!@name Iterators
//!@{
  using Iterator = std::vector<unsigned>::iterator;
  using ConstIterator = std::vector<unsigned>::const_iterator;
  inline Iterator begin() { return std::begin(map); }
  inline Iterator end() { return std::end(map); }
  inline ConstIterator begin() const { return std::begin(map); }
  inline ConstIterator end() const { return std::end(map); }
//!@}

  std::vector<unsigned> map;
};

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
 *   bond orders using Utils' BondDetector class.
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

/**
 * @brief Interpret molecules of a periodic system
 *
 * Interpret molecules into a simplified representation of structures with
 * periodic boundary conditions.
 *
 * @param atoms Atom collection with ghost atoms
 * @param periodicBonds Bond orders including extra bonds to ghost atoms
 * @param uninterestingAtoms List of atoms for which to skip shape
 *   classification and stereopermutator instantiation.
 * @param ghostAtomMap Map from ghost atom indices to their base atom
 * @param discretization Decide how bond orders are discretized into bond types
 * @param stereopermutatorThreshold If specified, limits the
 *   instantiation of BondStereopermutators onto edges whose fractional bond orders
 *   exceed the provided threshold. If this is not desired, specify boost::none.
 *
 * @pre Ghost atom indices must be a separate group of higher indices than both
 *   interesting and uninteresting atoms.
 *
 * @throws std::invalid_argument If the number of particles in the atom
 *   collection and bond order collection do not match.
 *
 * @parblock @note Assumes that the provided atom collection's positions are in
 * Bohr units @endparblock
 *
 * @parblock @note Any molecules interpreted with uninteresting atoms cannot be
 * passed to conformer generation routines. @endparblock
 *
 * @parblock @warning Ranking across periodic boundaries is incorrect.
 * @endparblock
 */
MASM_EXPORT MoleculesResult molecules(
  const Utils::AtomCollection& atoms,
  const Utils::BondOrderCollection& periodicBonds,
  const std::unordered_set<unsigned>& uninterestingAtoms,
  const std::unordered_map<unsigned, unsigned>& ghostAtomMap,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereopermutatorThreshold = 1.4
);

/*!
 * @brief Interpret molecules from Cartesian coordinates
 *
 * Attempts to interpret (possibly multiple) Molecules from element types and
 * positional information. Bond orders are calculated from atom-pairwise
 * spatial distances by van-der-Waals radii. The bond orders are then discretized
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
 */
MASM_EXPORT MoleculesResult molecules(
  const Utils::AtomCollection& atomCollection,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereopermutatorThreshold = 1.4
);

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

//! @brief Fn result datatype
struct FalsePositive {
  unsigned i;
  unsigned j;
  double probability;

  inline bool operator < (const FalsePositive& other) const {
    return probability < other.probability;
  }
};

/*! @brief Suggests false positives from a binary interpretation of bond orders
 *
 * Collects instances where shape classification is uncertain. If two connected
 * atoms are both uncertain (>= 50% probability that the point cloud could be
 * part of a random sample), lists the bond as a possible false postive.
 *
 * @warning Do not alter both bonds if there is a bond pair that have
 * overlapping indices. If suggested bonds overlap, remove only that bond with
 * the higher probabiltiy.
 */
MASM_EXPORT std::vector<FalsePositive> uncertainBonds(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders
);

/*! @brief Suggests false positive haptic ligand bonds from a binary interpretation of bond orders
 *
 * Generates a plane of best fit for each haptic ligand in the interpreted
 * graphs. If the angle of the normal of this plane to the axis defined by the
 * central atom and the site centroid is more than 30 degrees, tries to name a
 * single bond whose removal improves the interpretation.
 *
 * @note Suggested bonds can disconnect haptic sites. When making changes to
 * a bond order matrix based on suggestions from this function, apply them one
 * at a time based on the highest probability received. Additionally, if
 * multiple bonds must be removed to make a haptic ligand geometrically
 * reasonable, you will need to iteratively call this function and alter
 * suggested bond orders.
 */
MASM_EXPORT std::vector<FalsePositive> badHapticLigandBonds(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders
);

/*! @brief Iteratively applies false positive detection schemes
 *
 * @warning Pretty darn conservative implementation, removes only a single bond
 * from each false positive detection function call each iteration.
 */
MASM_EXPORT Utils::BondOrderCollection removeFalsePositives(
  const Utils::AtomCollection& atoms,
  Utils::BondOrderCollection bonds
);

} // namespace Interpret
} // namespace Molassembler
} // namespace Scine

#endif
