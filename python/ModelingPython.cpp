/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Modeling/BondDistance.h"

void init_modeling(pybind11::module& m) {
  using namespace Scine::Molassembler;

  auto modeling = m.def_submodule("modeling");
  modeling.doc() = "Modeling submodule";

  modeling.def(
    "bond_distance",
    &Bond::calculateBondDistance,
    pybind11::arg("a"),
    pybind11::arg("b"),
    pybind11::arg("bond_type"),
    "Calculates bond distance as modeled by UFF"
  );

  modeling.def(
    "bond_order",
    &Bond::calculateBondOrder,
    pybind11::arg("a"),
    pybind11::arg("b"),
    pybind11::arg("distance"),
    "Calculates bond order as modeled by UFF"
  );
}
