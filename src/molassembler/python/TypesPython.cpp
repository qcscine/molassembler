/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include "pybind11/stl.h"

#include "molassembler/Types.h"

void init_types(pybind11::module& m) {
  using namespace std::string_literals;
  using namespace Scine::molassembler;

  pybind11::enum_<BondType> bondType(
    m,
    "BondType",
    "Bond type enumeration. Besides the class organic single, double and triple "
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

  pybind11::class_<BondIndex> bondIndex(m, "BondIndex", "Ordered atom index pair");
  bondIndex.def(
    pybind11::init<AtomIndex, AtomIndex>(),
    pybind11::arg("a"),
    pybind11::arg("b"),
    "Initialize a bond index from two atom indices"
  );
  bondIndex.def_readwrite("first", &BondIndex::first);
  bondIndex.def_readwrite("second", &BondIndex::second);
  bondIndex.def(
    "__getitem__",
    [](const BondIndex& bond, unsigned index) -> AtomIndex {
      if(index > 1) {
        throw pybind11::index_error();
      }

      if(index == 0) {
        return bond.first;
      }

      return bond.second;
    }
  );
  bondIndex.def(
    "__setitem__",
    [](BondIndex& bond, unsigned index, AtomIndex i) {
      if(index > 1) {
        throw pybind11::index_error();
      }

      if(index == 0) {
        bond.first = i;
      } else {
        bond.second = i;
      }

      if(bond.first > bond.second) {
        std::swap(bond.first, bond.second);
      }
    }
  );
  bondIndex.def(
    "__contains__",
    [](const BondIndex& bond, AtomIndex i) {
      return (bond.first == i || bond.second == i);
    }
  );
  bondIndex.def(
    "__iter__",
    [](const BondIndex& bond) {
      return pybind11::make_iterator(
        bond.begin(),
        bond.end()
      );
    }
  );
  bondIndex.def(
    "__repr__",
    [](const BondIndex& bond) -> std::string {
      return "("s + std::to_string(bond.first) + ", "s + std::to_string(bond.second) + ")"s;
    }
  );
  bondIndex.def(pybind11::self == pybind11::self);
  bondIndex.def(pybind11::self < pybind11::self);

  /* AtomEnvironmentComponents binding
   * Cannot use None as an enum value since None is a reserved keyword in Python
   */
  pybind11::enum_<AtomEnvironmentComponents>(
    m,
    "AtomEnvironmentComponents",
    pybind11::arithmetic()
  ).value("NoComponents", AtomEnvironmentComponents::None)
   .value("ElementTypes", AtomEnvironmentComponents::ElementTypes)
   .value("BondOrders", AtomEnvironmentComponents::BondOrders)
   .value("Symmetries", AtomEnvironmentComponents::Symmetries)
   .value("Stereopermutations", AtomEnvironmentComponents::Stereopermutations)
   .value(
      "ElementsAndBonds",
      AtomEnvironmentComponents::ElementTypes
      | AtomEnvironmentComponents::BondOrders
    )
   .value(
      "ElementsBondsAndSymmetries",
      AtomEnvironmentComponents::ElementTypes
      | AtomEnvironmentComponents::BondOrders
      | AtomEnvironmentComponents::Symmetries
   )
   .value("All", AtomEnvironmentComponents::All);
}
