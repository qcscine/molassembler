/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Editing.h"
#include "Molassembler/Molecule.h"

void init_editing(pybind11::module& m) {
  using namespace Scine::Molassembler;

  auto editing = m.def_submodule("editing");
  editing.doc() = R"delim(
    A collection of functions to ease larger-scale molecule editing, since
    it can be difficult to get anywhere with the miniscule alteration functions
    defined in the :class:`Molecule` interface.
  )delim";

  editing.def(
    "cleave",
    pybind11::overload_cast<const Molecule&, BondIndex>(&Editing::cleave),
    pybind11::arg("molecule"),
    pybind11::arg("bridge"),
    R"delim(
      Cleave a molecule in two along a bridge bond.

      Bridge bonds are edges in the graph that whose removal splits the graph
      into two connected components. Any bonds in a cycle, for instance, are
      not bridge bonds.

      :param molecule: Molecule to cleave
      :param bridge: Bond index of bridge bond to cleave.
      :return: A pair of molecules
      :example:

      >>> import scine_utilities as utils
      >>> a = Molecule() # Makes H2
      >>> bond_index = a.addAtom(0, utils.ElementType.H) # Make linear H3
      >>> cleaved = editing.cleave(a, bond_index) # Split into H2 and H
    )delim"
  );

  editing.def(
    "cleave",
    pybind11::overload_cast<const Molecule&, Editing::AtomSitePair>(&Editing::cleave),
    pybind11::arg("molecule"),
    pybind11::arg("haptic_site"),
    R"delim(
      Cleave a molecule in two along a haptic site.

      Bridge bonds are edges in the graph that whose removal splits the graph
      into two connected components. Any bonds in a cycle, for instance, are
      not bridge bonds.

      :param molecule: Molecule to cleave
      :param haptic_site: Atom and site index pair indicating the haptic site
        to cleave
      :return: A pair of molecules. The first always contains the atom
        indicated by ``haptic_site``.
    )delim"
  );

  editing.def(
    "insert",
    &Editing::insert,
    pybind11::arg("log"),
    pybind11::arg("wedge"),
    pybind11::arg("log_bond"),
    pybind11::arg("first_wedge_atom"),
    pybind11::arg("second_wedge_atom"),
    R"delim(
      Insert a molecule into a bond of another molecule. Splits ``log`` at
      ``log_bond``, then inserts ``wedge`` at the split atoms, connecting the
      first atom of ``log_bond`` with ``first_wedge_atom`` and the second
      ``log_bond`` atom with ``second_wedge_atom``.

      The bond type of the ``log_bond`` is reused in the new bonds formed to the
      ``wedge`` atoms.

      :param log: The molecule being inserted into
      :param wedge: The molecule being inserted into the ``log``
      :param log_bond: Log's bond that ``wedge`` should be inserted into
      :param first_wedge_atom: The atom of ``wedge`` to bond to the first atom
        in ``log_bond``
      :param second_wedge_atom: The atom of ``wedge`` to bond to the second
        atom in ``log_bond``
      :return: The result of the insert operation
    )delim"
  );

  editing.def(
    "superpose",
    &Editing::superpose,
    pybind11::arg("top"),
    pybind11::arg("bottom"),
    pybind11::arg("top_overlay_atom"),
    pybind11::arg("bottom_overlay_atom"),
    R"delim(
      Fuse two molecules, adding all adjacencies of one Molecule's atoms to
      another

      Adds all adjacent atoms and continuations of ``bottom_overlay_atom`` in
      bottom to ``top_overlay_atom`` in top. ``top_overlay_atom``'s element
      type is unchanged as it is the 'top' of the superimposition / overlay.

      :param top: The molecule at the "top" of the superposition.
      :param bottom: The molecule at the "bottom" of the superposition.
      :param top_overlay_atom: The atom of ``top`` that is placed "onto"
        ``bottom``'s ``bottom_overlay_atom``
      :param bottom_overlay_atom: The atom of ``bottom`` to place "beneath"
        top's ``top_overlay_atom``
    )delim"
  );

  editing.def(
    "substitute",
    &Editing::substitute,
    pybind11::arg("left"),
    pybind11::arg("right"),
    pybind11::arg("left_bridge"),
    pybind11::arg("right_bridge"),
    R"delim(
      Connect two molecules by substituting away the lighter side of a pair of
      bonds of separate molecules.

      The heavy side is chosen by number of atoms first, then molecular weight
      if the number of atoms is equal. Should both sides be equal in both, which
      side is picked is undefined.

      :param left: The first molecule
      :param right: The second molecule
      :param left_bridge: Left's bridge bond from which to substitute the
        lighter part away.
      :param right_bridge: Right's bridge bond from which to substitute the
        lighter part away.
    )delim"
  );

  editing.def(
    "connect",
    &Editing::connect,
    pybind11::arg("left"),
    pybind11::arg("right"),
    pybind11::arg("left_atom"),
    pybind11::arg("right_atom"),
    pybind11::arg("bond_type"),
    R"delim(
      Connect two molecules by creating a new bond between two atoms from
      separate molecules

      :param left: The first molecule
      :param right: The second molecule
      :param left_atom: The atom from ``left`` to connect
      :param right_atom: The atom from ``right`` to connect
      :param bond_type: The bond type with which to connect ``left_atom`` and
        ``right_atom``
    )delim"
  );

  editing.def(
    "add_ligand",
    &Editing::addLigand,
    pybind11::arg("a"),
    pybind11::arg("ligand"),
    pybind11::arg("complexating_atom"),
    pybind11::arg("ligand_binding_atoms"),
    R"delim(
      Connect two molecules by connecting multiple atoms from one to a single
      atom of the other via single bonds.

      :param a: The molecule the ligand is being connected to
      :param ligand: The ligand molecule being bound
      :param complexating_atom: The atom in ``a`` to bind ligand to
      :param ligand_binding_atoms: Atoms in ``ligand`` to bind to
        ``complexating_atom`` to.
    )delim"
  );
}
