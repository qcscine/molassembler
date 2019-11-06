/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/operators.h"
#include "OptionalPython.h"

#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/StereopermutatorList.h"

#include "Utils/Geometry/ElementTypes.h"
#include "Utils/Geometry/AtomCollection.h"

#include "boost/process/child.hpp"
#include "boost/process/io.hpp"
#include "boost/process/search_path.hpp"

bool graphvizInPath() {
  return !boost::process::search_path("dot").empty();
}

std::string pipeSVG(const Scine::molassembler::Molecule& molecule) {
  std::string callString = "dot -Tsvg";
  std::stringstream os;

  // Construct pipe streams for redirection
  boost::process::opstream ips;
  boost::process::pstream ps;
  boost::process::pstream err;

  // Start the child process
  boost::process::child childProcess(callString, boost::process::std_in<ips, boost::process::std_out> ps,
                                     boost::process::std_err > err);

  // Feed our graphviz into the process
  ips << molecule.dumpGraphviz();
  ips.flush();
  ips.pipe().close();

  // Wait for the child process to exit
  childProcess.wait();

  std::stringstream stderrStream;
#if BOOST_VERSION >= 107000
  /* NOTE: This implementation of buffer transfers in boost process has a bug
   * that isn't fixed before Boost 1.70.
   */
  os << ps.rdbuf();
  stderrStream << err.rdbuf();
#else
  // Workaround: cast to a parent class implementing rdbuf() correctly.
  using BasicIOSReference = std::basic_ios<char, std::char_traits<char>>&;
  // Feed the results into our ostream
  os << static_cast<BasicIOSReference>(ps).rdbuf();
  stderrStream << static_cast<BasicIOSReference>(err).rdbuf();
#endif

  return os.str();
}

void init_molecule(pybind11::module& m) {
  using namespace Scine::molassembler;
  using namespace Scine::Utils;

  pybind11::class_<Molecule> molecule(
    m,
    "Molecule",
    "Models a molecule as a graph and a list of stereopermutators"
  );

  /* Constructors */
  molecule.def(
    pybind11::init<>(),
    "Initialize a hydrogen molecule"
  );

  molecule.def(
    pybind11::init<ElementType>(),
    "Initialize a single-atom molecule"
  );

  molecule.def(
    pybind11::init<ElementType, ElementType, BondType>(),
    pybind11::arg("first_element"),
    pybind11::arg("second_element"),
    pybind11::arg("bond_type") = BondType::Single,
    "Initialize a molecule from two element types and a mutual bond type"
  );

  molecule.def(
    pybind11::init<OuterGraph>(),
    pybind11::arg("graph"),
    R"delim(
      Initialize a molecule from connectivity alone, inferring shapes and
      stereopermutators from the graph.
    )delim"
  );

  molecule.def(
    "hash",
    &Molecule::hash,
    R"delim(
      Calculates a convoluted hash of a molecule. The molecule must be at least
      partially canonical. Hashes between molecules of different canonicity are
      not comparable.
    )delim"
  );

  molecule.def(
    "__hash__",
    &Molecule::hash
  );

  molecule.def_static(
    "apply_canonicalization_map",
    &Molecule::applyCanonicalizationMap,
    pybind11::arg("canonicalization_index_map"),
    pybind11::arg("atom_collection"),
    R"delim(
      Reorders an atom collection according to an index mapping from
      canonicalization

      :param canonicalization_index_map: Index mapping saved from previous
        canonicalization
      :param atom_collection: Atom collection to reorder
      :return: Reordered atom collection
    )delim"
  );

  /* Modifiers */
  molecule.def(
    "add_atom",
    &Molecule::addAtom,
    pybind11::arg("element"),
    pybind11::arg("adjacent_to"),
    pybind11::arg("bond_type") = BondType::Single,
    R"delim(
      Add an atom to the molecule, attaching it to an existing atom by a
      specified bond type.

      :param element: Element type of the new atom
      :param adjacent_to: Atom to which the new atom is added
      :param bond_type: Bond type with which the new atom is attached
    )delim"
  );

  molecule.def(
    "add_bond",
    &Molecule::addBond,
    pybind11::arg("first_atom"),
    pybind11::arg("second_atom"),
    pybind11::arg("bond_type") = BondType::Single,
    R"delim(
      Adds a bond between two existing atoms.

      :param first_atom: First atom to bond
      :param second_atom: Second atom to bond
      :param bond_type: Bond type with which to bond the atoms
    )delim"
  );

  molecule.def(
    "assign_stereopermutator",
    pybind11::overload_cast<AtomIndex, const boost::optional<unsigned>&>(
      &Molecule::assignStereopermutator
    ),
    pybind11::arg("atom"),
    pybind11::arg("assignment_option"),
    R"delim(
      Sets the atom stereopermutator assignment at a particular atom

      :param atom: Atom index of the atom stereopermutator to set
      :param assignment_option: An assignment integer if the stereopermutator
        is to be assigned or None if the stereopermutator is to be dis-assigned.
    )delim"
  );

  molecule.def(
    "assign_stereopermutator",
    pybind11::overload_cast<const BondIndex&, const boost::optional<unsigned>&>(
      &Molecule::assignStereopermutator
    ),
    pybind11::arg("bond_index"),
    pybind11::arg("assignment_option"),
    R"delim(
      Sets the bond stereopermutator assignment at a particular bond

      :param bond_index: Bond index of the stereopermutator to set
      :param assignment_option: An assignment integer if the stereopermutator
        is to be assigned or None if the stereopermutator is to be dis-assigned.
    )delim"
  );

  molecule.def(
    "assign_stereopermutator_randomly",
    [](Molecule& mol, AtomIndex a) {
      mol.assignStereopermutatorRandomly(a, randomnessEngine());
    },
    pybind11::arg("atom"),
    R"delim(
      Assigns an atom stereopermutator at random (assignments are weighted by
      relative statistical occurence).

      :param atom: Atom index of the stereopermutator to assign randomly.
    )delim"
  );

  molecule.def(
    "assign_stereopermutator_randomly",
    [](Molecule& mol, const BondIndex& b) {
      mol.assignStereopermutatorRandomly(b, randomnessEngine());
    },
    pybind11::arg("bond_index"),
    R"delim(
      Assigns a bond stereopermutator at random.

      :param bond_index: Bond index of the stereopermutator to assign randomly.
    )delim"
  );

  molecule.def(
    "canonicalize",
    &Molecule::canonicalize,
    pybind11::arg("components_bitmask") = AtomEnvironmentComponents::All,
    R"delim(
      Transform the molecule to a canonical form. Invalidates all atom and bond
      indices.

      :param components_bitmask: The components of the molecular graph to
        include in the canonicalization procedure.
      :return: Flat index mapping/permutation from old indices to new
    )delim"
  );

  molecule.def(
    "remove_atom",
    &Molecule::removeAtom,
    pybind11::arg("atom"),
    R"delim(
      Remove an atom from the graph, including bonds to it, after checking
      that removing it is safe, i.e. the removal does not disconnect the graph.

      :param atom: Atom to remove
    )delim"
  );

  molecule.def(
    "remove_bond",
    pybind11::overload_cast<AtomIndex, AtomIndex>(&Molecule::removeBond),
    pybind11::arg("first_atom"),
    pybind11::arg("second_atom"),
    R"delim(
      Remove a bond from the graph, after checking that removing it is safe,
      i.e. the removal does not disconnect the graph.

      :param first_atom: First atom of the bond to be removed
      :param second_atom: Second atom of the bond to be removed
    )delim"
  );

  molecule.def(
    "remove_bond",
    pybind11::overload_cast<const BondIndex&>(&Molecule::removeBond),
    pybind11::arg("bond_index"),
    R"delim(
      Remove a bond from the graph, after checking that removing it is safe,
      i.e. the removal does not disconnect the graph.

      :param bond_index: Bond index of the bond to be removed
    )delim"
  );

  molecule.def(
    "set_bond_type",
    &Molecule::setBondType,
    pybind11::arg("first_atom"),
    pybind11::arg("second_atom"),
    pybind11::arg("bond_type"),
    R"delim(
      Change the bond type between two atoms. Inserts the bond if it doesn't
      yet exist.

      :param first_atom: First atom of the bond to be changed
      :param second_atom: Second atom of the bond to be changed
      :param bond_type: The new bond type
      :return: Whether the bond already existed
    )delim"
  );

  molecule.def(
    "set_element_type",
    &Molecule::setElementType,
    pybind11::arg("atom"),
    pybind11::arg("element"),
    R"delim(
      Change the element type of an atom.

      :param atom: Atom index of the atom to alter
      :param element: New element type to set
    )delim"
  );

  molecule.def(
    "set_shape_at_atom",
    &Molecule::setShapeAtAtom,
    pybind11::arg("atom"),
    pybind11::arg("shape"),
    R"delim(
      Change the local shape at an atom.

      This sets the local shape at a specific atom index. There are a number of
      cases that this function treats differently, besides faulty arguments: If
      there is already a AtomStereopermutator instantiated at this atom index,
      its underlying shape is altered. If there is no AtomStereopermutator at
      this index, one is instantiated. In all cases, new or modified
      stereopermutators are default-assigned if there is only one possible
      assignment.
    )delim"
  );

  /* Information */
  molecule.def(
    "dump_graphviz",
    &Molecule::dumpGraphviz,
    "Returns a graphviz string representation of the molecule"
  );

  molecule.def_property_readonly(
    "graph",
    &Molecule::graph,
    "Read only access to the graph representation"
  );

  molecule.def_property_readonly(
    "stereopermutators",
    &Molecule::stereopermutators,
    "Read only access to the list of stereopermutators"
  );

  molecule.def_property_readonly(
    "canonical_components",
    &Molecule::canonicalComponents,
    "Yields the components of the molecule that were used in a previous canonicalization. Can be None."
  );

  molecule.def(
    "canonical_compare",
    &Molecule::canonicalCompare,
    pybind11::arg("other"),
    pybind11::arg("components_bitmask") = AtomEnvironmentComponents::All,
    R"delim(
      Modular comparison of this Molecule with another, assuming that both are
      in a canonical form.

      :param other: The other (canonical) molecule to compare against
      :param components_bitmask: The components of an atom's environment to
        include in the comparison. You should use the same bitmask as when
        canonicalizing the molecules you are comparing here. It may be possible
        to use a bitmask with fewer components, but certainly not one with more.
    )delim"
  );

  molecule.def(
    "partial_compare",
    &Molecule::modularCompare,
    pybind11::arg("other"),
    pybind11::arg("components_bitmask"),
    R"delim(
      Modular comparison of this Molecule with another.

      This permits detailed specification of which elements of the molecular
      information you want to use in the comparison.

      Equality comparison is performed in several stages: First, at each atom
      position, a hash is computed that encompasses all local information that
      is specified to be used by the components_bitmask parameter. This hash is
      then used during graph isomorphism calculation to avoid finding an
      isomorphism that does not consider the specified factors.

      If an isomorphism is found, it is then validated. Bond orders and
      stereopermutators across both molecules are compared using the found
      isomorphism as an index map.

      Note that this function is not faster for molecules stored in any
      (possibly partially) canonical form. Use canonical_compare for molecules
      that have been canonicalized to some degree. Note also that equality
      comparison defaults to fast comparisons if both instances are fully
      canonical.

      :param other: The molecule to compare against
      :param components_bitmask: The components of the molecule to use in the
        comparison
    )delim"
  );

  molecule.def(pybind11::self == pybind11::self);
  molecule.def(pybind11::self != pybind11::self);

  /* Integration with IPython / Jupyter */
  molecule.def(
    "_repr_svg_",
    &::pipeSVG,
    "Generates an SVG representation of the molecule"
  );
}
