/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

void init_atom_stereopermutator(pybind11::module& m);
void init_bond_stereopermutator(pybind11::module& m);
void init_conformers(pybind11::module& m);
void init_cycles(pybind11::module& m);
void init_directed_conformer_generator(pybind11::module& m);
void init_editing(pybind11::module& m);
void init_interpret(pybind11::module& m);
void init_io(pybind11::module& m);
void init_molecule(pybind11::module& m);
void init_options(pybind11::module& m);
void init_graph(pybind11::module& m);
void init_random_engine(pybind11::module& m);
void init_ranking_information(pybind11::module& m);
void init_serialization(pybind11::module& m);
void init_stereopermutator_list(pybind11::module& m);
void init_shape_submodule(pybind11::module& m);
void init_types(pybind11::module& m);
void init_version(pybind11::module& m);

PYBIND11_MODULE(scine_molassembler, m) {
  m.doc() = R"(
    Pybind11 bindings for scine_molassembler

    .. currentmodule:: scine_molassembler

    .. autosummary::
       :toctree:
  )";

  // Requires other modules to function properly
  // auto utils = pybind11::module::import("scine_utilities");

  // Order is important here, do not reorder
  init_version(m);
  init_types(m);
  init_random_engine(m);
  init_options(m);
  init_shape_submodule(m);
  init_cycles(m);
  init_graph(m);
  init_ranking_information(m);
  init_atom_stereopermutator(m);
  init_bond_stereopermutator(m);
  init_stereopermutator_list(m);
  init_molecule(m);
  init_editing(m);
  init_interpret(m);
  init_io(m);
  init_serialization(m);
  init_conformers(m);
  init_directed_conformer_generator(m);
}
