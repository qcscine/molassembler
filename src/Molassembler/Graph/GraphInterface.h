/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Interface class for the molecular graph
 */

#ifndef INCLUDE_MOLASSEMBLER_GRAPH_INTERFACE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_INTERFACE_H

#include "Molassembler/Types.h"

#include "boost/optional/optional_fwd.hpp"
#include "Utils/Geometry/ElementTypes.h"

namespace Scine {
namespace Molassembler {

class Cycles;

struct GraphInterface {
  virtual ~GraphInterface() = default;

//!@name Modification
//!@{
  virtual AtomIndex addAtom(Utils::ElementType e, AtomIndex i, BondType type) = 0;

  virtual BondIndex addBond(AtomIndex i, AtomIndex j, BondType type) = 0;

  virtual void setElementType(AtomIndex i, Utils::ElementType type) = 0;

  virtual bool setBondType(AtomIndex i, AtomIndex j, BondType type) = 0;

  virtual void removeAtom(AtomIndex i) = 0;

  virtual void removeBond(const BondIndex& bond) = 0;
//!@}

//!@name Information
//!@{
  virtual bool adjacent(AtomIndex a, AtomIndex b) const = 0;

  virtual boost::optional<BondIndex> bond(AtomIndex a, AtomIndex b) const = 0;

  virtual BondType bondType(const BondIndex& edge) const = 0;

  virtual bool canRemove(AtomIndex a) const = 0;

  virtual bool canRemove(const BondIndex& edge) const = 0;

  virtual const Cycles& cycles() const = 0;

  virtual unsigned degree(AtomIndex a) const = 0;

  virtual Utils::ElementType elementType(AtomIndex a) const = 0;

  virtual AtomIndex V() const = 0;

  virtual unsigned E() const = 0;
//!@}
};

} // namespace Molassembler
} // namespace Scine

#endif
