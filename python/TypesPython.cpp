/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/operators.h"

#include "Molassembler/Types.h"

void init_types(pybind11::module& m) {
  using namespace std::string_literals;
  using namespace Scine::Molassembler;

  pybind11::enum_<BondType> bondType(
    m,
    "BondType",
    R"delim(
      Bond type enumeration. Besides the classic organic single, double and
      triple bonds, bond orders up to sextuple are explicitly included. Eta is
      a bond order used internally by the library to represent haptic bonding.
      It should not be set by users.
    )delim"
  );

  bondType
  .value("Single", BondType::Single, "Single bond")
  .value("Double", BondType::Double, "Double bond")
  .value("Triple", BondType::Triple, "Triple bond")
  .value("Quadruple", BondType::Quadruple, "Quadruple bond")
  .value("Quintuple", BondType::Quintuple, "Quintuple bond")
  .value("Sextuple", BondType::Sextuple, "Sextuple bond")
  .value("Eta", BondType::Eta, "Eta bond, indicates haptic bonding");

  // Leave out LengthUnit, it should not be necessary

  pybind11::class_<BondIndex> bondIndex(
    m,
    "BondIndex",
    R"delim(
      Ordered atom index pair indicating a bond in a molecule graph.

      Has some container magic methods and comparators.

      :example:

      >>> b = BondIndex(10, 4)
      >>> b
      (4, 10)
      >>> b.first
      4
      >>> b[1]
      10
      >>> b[1] = 3 # Change part of the bond index
      >>> b # The altered object is still ordered
      (3, 4)
      >>> 3 in b
      True
      >>> 10 in b
      False
      >>> c = BondIndex(4, 3)
      >>> b == c
      True
      >>> d = BondIndex(2, 4)
      >>> d < c
      True
    )delim"
  );
  bondIndex.def(
    pybind11::init<AtomIndex, AtomIndex>(),
    pybind11::arg("a"),
    pybind11::arg("b"),
    "Initialize a bond index from two atom indices"
  );
  bondIndex.def_readwrite("first", &BondIndex::first, "The lesser atom index");
  bondIndex.def_readwrite("second", &BondIndex::second, "The greater atom index");
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
    pybind11::arithmetic(),
    R"delim(
      Denotes information parts of molecules. Relevant for molecular comparison
      and hashing.
    )delim"
  ).value("Connectivity", AtomEnvironmentComponents::Connectivity, "Consider only the graph")
   .value("ElementTypes", AtomEnvironmentComponents::ElementTypes, "Element types")
   .value("BondOrders", AtomEnvironmentComponents::BondOrders, "Bond orders")
   .value("Shapes", AtomEnvironmentComponents::Shapes, "Shapes")
   .value("Stereopermutations", AtomEnvironmentComponents::Stereopermutations, "Stereopermutations")
   .value(
      "ElementsAndBonds",
      AtomEnvironmentComponents::ElementTypes
      | AtomEnvironmentComponents::BondOrders,
      "Consider element types and bond orders"
    )
   .value(
      "ElementsBondsAndShapes",
      AtomEnvironmentComponents::ElementTypes
      | AtomEnvironmentComponents::BondOrders
      | AtomEnvironmentComponents::Shapes,
      "Consider element types, bond orders and shapes"
   )
   .value("All", AtomEnvironmentComponents::All, "Consider element types, bond orders, shapes and stereopermutations");
}
