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
    "Cleave a molecule in two along a bridge bond"
  );

  editing.def(
    "insert",
    &Editing::insert,
    "Insert a molecule into a bond of another molecule"
  );

  editing.def(
    "superpose",
    &Editing::superpose,
    "Fuse two molecules, adding all adjacencies of one Molecule's atoms to "
    "another"
  );

  editing.def(
    "substitute",
    &Editing::substitute,
    "Connect two molecules by substituting away the lighter side of a pair of "
    "bonds of separate molecules"
  );

  editing.def(
    "connect",
    &Editing::connect,
    "Connect two molecules by creating a new bond between two atoms from "
    "separate molecules"
  );
}
