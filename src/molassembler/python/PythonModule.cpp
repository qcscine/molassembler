/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"

void init_atom_stereopermutator(pybind11::module& m);
void init_bond_stereopermutator(pybind11::module& m);
void init_conformers(pybind11::module& m);
void init_cycles(pybind11::module& m);
void init_interpret(pybind11::module& m);
void init_io(pybind11::module& m);
void init_molecule(pybind11::module& m);
void init_options(pybind11::module& m);
void init_outer_graph(pybind11::module& m);
void init_random_engine(pybind11::module& m);
void init_ranking_information(pybind11::module& m);
void init_serialization(pybind11::module& m);
void init_stereopermutator_list(pybind11::module& m);
void init_symmetry_submodule(pybind11::module& m);
void init_types(pybind11::module& m);

PYBIND11_MODULE(molassembler, m) {
  m.doc() = "Pybind11 Bindings for molassembler";

  // Order is important here, do not reorder
  init_types(m);
  init_random_engine(m);
  init_options(m);
  init_symmetry_submodule(m);
  init_cycles(m);
  init_outer_graph(m);
  init_ranking_information(m);
  init_atom_stereopermutator(m);
  init_bond_stereopermutator(m);
  init_stereopermutator_list(m);
  init_molecule(m);
  init_interpret(m);
  init_io(m);
  init_serialization(m);
  init_conformers(m);
}
