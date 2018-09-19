// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_INTERPRET_H
#define INCLUDE_MOLASSEMBLER_INTERPRET_H

#include "boost/optional.hpp"
#include <vector>

/*!@file
 *
 * @brief Interpret multiple molecules in positional information
 *
 * Contains functionality permitting the interpretation of several Molecules
 * from three-dimensional structures with or without accompanying bond orders.
 */

// External forward-declarations
namespace Delib {
class AtomCollection;
class BondOrderCollection;
class ElementTypeCollection;
} // namespace Delib

namespace molassembler {

// Forward-declarations
class Molecule;
class AngstromWrapper;

//! Specify the algorithm used to discretize floating-point bond orders into bond types
enum class BondDiscretizationOption : unsigned {
  //! All bond orders >= 0.5 are considered single bonds
  Binary,
  //! Bond orders are rounded to the nearest integer
  RoundToNearest
};

//! Result type of an interpret call.
struct InterpretResult {
  //! The list of individual molecules found in the 3D information
  std::vector<Molecule> molecules;
  //! A map from an index within the AtomCollection to its index in molecules member
  std::vector<unsigned> componentMap;
};

//! The function that actually does all the work with the library-internal wrapper
InterpretResult interpret(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const Delib::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereocenterBondOrderThresholdOptional = 1.4
);

//! A forwarding function with the library-internal wrapper types
InterpretResult interpret(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereocenterBondOrderThresholdOptional = 1.4
);

/*! Interpret molecules from element types, positional information and a bond order collection.
 *
 * Attempts to interpret (possibly multiple) Molecules from element types,
 * positional information and a bond order collection. Bond orders are
 * discretized into bond types. Connected components within the space are
 * identified and individually instantiated into Molecules. The instantiation
 * behavior of BondStereocenters in the Molecules can be limited to edges whose
 * bond order exceeds a particular value.
 *
 * \param atomCollection Element types and positional information in Bohr units.
 * \param bondOrders Fractional bond orders
 * \param discretization Decide how bond orders are discretized into bond types
 * \param stereocenterBondOrderThresholdOptional If specified, limits the
 *   instantiation of BondStereocenters onto edges whose fractional bond orders
 *   exceed the provided threshold. If this is not desired, specify boost::none.
 *
 * \note Assumes that the provided atom collection's positions are in
 * Bohr units.
 */
InterpretResult interpret(
  const Delib::AtomCollection& atomCollection,
  const Delib::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereocenterBondOrderThresholdOptional = 1.4
);

/*! Interpret molecules in 3D information.
 *
 * Attempts to interpret (possibly multiple) Molecules from element types and
 * positional information. Bond orders are calculated from atom-pairwise
 * spatial distances using UFF parameters. The bond orders are then discretized
 * into bond types. Connected components within the space are identified and
 * individually instantiated into Molecules. The instantiation behavior of
 * BondStereocenters in the Molecules can be limited to edges whose bond order
 * exceeds a particular value.
 *
 * \param atomCollection Element types and positional information in Bohr units.
 * \param bondOrders Fractional bond orders
 * \param discretization Decide how bond orders are discretized into bond types
 * \param stereocenterBondOrderThresholdOptional If specified, limits the
 *   instantiation of BondStereocenters onto edges whose fractional bond orders
 *   exceed the provided threshold
 *
 * \note Assumes that the provided atom collection's positions are in
 * Bohr units.
 *
 * \warning UFF parameter bond order calculation is very primitive and carries
 *   a high risk of misinterpretation
 */
InterpretResult interpret(
  const Delib::AtomCollection& atomCollection,
  BondDiscretizationOption discretization = BondDiscretizationOption::Binary,
  const boost::optional<double>& stereocenterBondOrderThresholdOptional = 1.4
);

} // namespace molassembler

#endif
