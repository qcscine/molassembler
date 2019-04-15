/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "molassembler/Editing.h"
#include "molassembler/Molecule.h"

void init_editing(pybind11::module& m) {
  using namespace Scine::molassembler;

  auto editing = m.def_submodule("editing");

  editing.def(
    "cleave",
    &Editing::cleave,
    pybind11::arg("molecule"),
    pybind11::arg("bridge"),
    "Cleave a molecule in two along a bridge bond"
  );

  editing.def(
    "insert",
    &Editing::insert,
    pybind11::arg("log"),
    pybind11::arg("wedge"),
    pybind11::arg("log_bond"),
    pybind11::arg("first_wedge_atom"),
    pybind11::arg("second_wedge_atom"),
    "Insert a molecule into a bond of another molecule. Splits log at log_bond, "
    "then inserts wedge at the split atoms, connecting the first atom of "
    "log_bond with first_wedge_atom and the second log_bond atom with "
    "second_wedge_atom."
  );

  editing.def(
    "superpose",
    &Editing::superpose,
    pybind11::arg("top"),
    pybind11::arg("bottom"),
    pybind11::arg("top_overlay_atom"),
    pybind11::arg("bottom_overlay_atom"),
    "Fuse two molecules, adding all adjacencies of one Molecule's atoms to "
    "another"
  );

  editing.def(
    "substitute",
    &Editing::substitute,
    pybind11::arg("left"),
    pybind11::arg("right"),
    pybind11::arg("right_bridge"),
    pybind11::arg("left_bridge"),
    "Connect two molecules by substituting away the lighter side of a pair of "
    "bonds of separate molecules"
  );

  editing.def(
    "connect",
    &Editing::connect,
    pybind11::arg("left"),
    pybind11::arg("right"),
    pybind11::arg("left_atom"),
    pybind11::arg("right_atom"),
    pybind11::arg("bond_type"),
    "Connect two molecules by creating a new bond between two atoms from "
    "separate molecules"
  );

  editing.def(
    "add_ligand",
    &Editing::addLigand,
    pybind11::arg("a"),
    pybind11::arg("ligand"),
    pybind11::arg("complexating_atom"),
    pybind11::arg("ligand_binding_atoms"),
    "Connect two molecules by connecting multiple atoms from one to a single "
    "atom of the other via single bonds"
  );
}
