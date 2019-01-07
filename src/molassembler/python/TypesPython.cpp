/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"

#include "molassembler/Types.h"

void init_types(pybind11::module& m) {
  using namespace Scine::molassembler;
  pybind11::enum_<BondType> bondType(
    m,
    "BondType",
    "Bond type numeration. Besides the class organic single, double and triple "
    "bonds, bond orders up to sextuple are explicitly included. Eta is a bond "
    "order used internally by the library to represent haptic bonding. It "
    "should not be set by users."
  );

  bondType
  .value("Single", BondType::Single)
  .value("Double", BondType::Double)
  .value("Triple", BondType::Triple)
  .value("Quadruple", BondType::Quadruple)
  .value("Quintuple", BondType::Quintuple)
  .value("Sextuple", BondType::Sextuple)
  .value("Eta", BondType::Eta);

  // Leave out LengthUnit, it should not be necessary

  pybind11::class_<BondIndex> bondIndex(m, "BondIndex", "Bond index with ordered atom indices");
  bondIndex.def_readwrite("first", &BondIndex::first);
  bondIndex.def_readwrite("second", &BondIndex::second);
  bondIndex.def(pybind11::self == pybind11::self);
  bondIndex.def(pybind11::self < pybind11::self);
}
