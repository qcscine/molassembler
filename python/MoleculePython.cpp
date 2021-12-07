/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/operators.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/Serialization.h"

#include "Utils/Geometry/ElementTypes.h"
#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Geometry/FormulaGenerator.h"

#include "boost/process/child.hpp"
#include "boost/process/io.hpp"
#include "boost/process/search_path.hpp"

namespace {

bool graphvizInPath() {
  return !boost::process::search_path("dot").empty();
}

std::string pipeSVG(const Scine::Molassembler::Molecule& molecule) {
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

} // namespace

void init_molecule(pybind11::module& m) {
  using namespace Scine::Molassembler;
  using namespace Scine::Utils;

  pybind11::class_<Molecule> molecule(
    m,
    "Molecule",
    R"delim(
      Models a molecule as a :class:`Graph` and a :class:`StereopermutatorList`.
    )delim"
  );

  /* Constructors */
  molecule.def(
    pybind11::init<>(),
    R"delim(
      Initialize a hydrogen molecule

      >>> h2 = Molecule()
      >>> h2.graph.V
      2
      >>> h2.graph.E
      1
    )delim"
  );

  molecule.def(
    pybind11::init<ElementType>(),
    R"delim(
      Initialize a single-atom molecule.

      This is a bit of a paradox, yes, and it might have been preferable for
      the concept of a molecule to contain at least two bonded atoms, but
      unfortunately single atoms occur everywhere and enforcing the concept
      would complicate many interfaces.

      >>> import scine_utilities as utils
      >>> f = Molecule(utils.ElementType.F)
      >>> f.graph.V
      1
      >>> f.graph.E
      0
    )delim"
  );

  molecule.def(
    pybind11::init<ElementType, ElementType, BondType>(),
    pybind11::arg("first_element"),
    pybind11::arg("second_element"),
    pybind11::arg("bond_type") = BondType::Single,
    R"delim(
      Initialize a molecule from two element types and a mutual :class:`BondType`

      >>> # Make H-F
      >>> import scine_utilities as utils
      >>> hf = Molecule(utils.ElementType.H, utils.ElementType.F)
      >>> hf.graph.V == 2
      True
    )delim"
  );

  molecule.def(
    pybind11::init<Graph>(),
    pybind11::arg("graph"),
    R"delim(
      Initialize a molecule from connectivity alone, inferring shapes and
      stereopermutators from the graph.

      >>> # Rebuild a molecule with an assigned stereopermutator from just the graph
      >>> a = io.experimental.from_smiles("[C@](F)(Cl)(C)[H]")
      >>> a.stereopermutators.has_unassigned_permutators()
      False
      >>> b = Molecule(a.graph)
      >>> b.stereopermutators.has_unassigned_permutators()
      True
    )delim"
  );

  molecule.def(
    "hash",
    &Molecule::hash,
    R"delim(
      Calculates a convoluted hash of a molecule. The molecule must be at least
      partially canonical. Hashes between molecules of different canonicity are
      not comparable.

      >>> # Show that hash values differ at various levels of canonicity
      >>> from copy import copy
      >>> spiro = io.experimental.from_smiles("C12(CCC1)CCC2")
      >>> # We make two variants of the molecule that have different canonicalization states
      >>> # to demonstrate that their hashes are unequal. We discard the mappings
      >>> # we get from canonicalize()
      >>> partially_canonical = copy(spiro)
      >>> _ = partially_canonical.canonicalize(AtomEnvironmentComponents.ElementsAndBonds)
      >>> fully_canonical = copy(spiro)
      >>> _ = fully_canonical.canonicalize()
      >>> partially_canonical == fully_canonical
      True
      >>> partially_canonical.hash() == fully_canonical.hash()
      False
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
      canonicalization.

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
      :param bond_type: :class:`BondType` with which the new atom is attached

      >>> # Let's make linear H3
      >>> import scine_utilities as utils
      >>> mol = Molecule() # Default constructor makes H2
      >>> _ = mol.add_atom(utils.ElementType.H, 0) # Make linear H3
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
      :param bond_type: :class:`BondType` with which to bond the atoms

      >>> # Let's make triangular H3
      >>> import scine_utilities as utils
      >>> mol = Molecule() # Default constructor makes H2
      >>> _ = mol.add_atom(utils.ElementType.H, 0) # Make linear H3
      >>> _ = mol.add_bond(1, 2) # Make triangular H3
    )delim"
  );

  molecule.def(
    "add_permutator",
    &Molecule::addPermutator,
    pybind11::arg("bond"),
    pybind11::arg("alignment") = BondStereopermutator::Alignment::Eclipsed,
    R"delim(
      Add a BondStereopermutator to the molecule

      .. note:: You can't add AtomStereopermutators to the molecule manually.
         These are automatically present on non-terminal atoms.

      :param bond: Bond to place the permutator at
      :param alignment: Alignment with which to construct the permutator

      :returns: A reference to the added stereopermutator
      :raises RuntimeError: If there is already a permutator present at this
        bond
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

      :param atom: Atom index of the :class:`AtomStereopermutator` to set
      :param assignment_option: An assignment integer if the stereopermutator
        is to be assigned or ``None`` if the stereopermutator is to be dis-assigned.

      >>> # Assign an unspecified asymmetric carbon atom and then dis-assign it
      >>> mol = io.experimental.from_smiles("F[CH1](Br)C")
      >>> asymmetric_carbon_index = 1
      >>> mol.assign_stereopermutator(asymmetric_carbon_index, 0)
      >>> mol.stereopermutators.option(asymmetric_carbon_index).assigned
      0
      >>> mol.assign_stereopermutator(asymmetric_carbon_index, None)
      >>> mol.stereopermutators.option(asymmetric_carbon_index).assigned is None
      True
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

      :param bond_index: :class:`BondIndex` of the :class:`BondStereopermutator` to set
      :param assignment_option: An assignment integer if the stereopermutator
        is to be assigned or ``None`` if the stereopermutator is to be
        dis-assigned.

      >>> # Dis-assign an assigned bond stereopermutator
      >>> ethene = io.experimental.from_smiles("C/C=C\C")
      >>> double_bond_index = BondIndex(1, 2)
      >>> assert ethene.graph.bond_type(double_bond_index) == BondType.Double
      >>> ethene.stereopermutators.option(double_bond_index).assigned is not None
      True
      >>> ethene.assign_stereopermutator(double_bond_index, None)
      >>> ethene.stereopermutators.option(double_bond_index).assigned is not None
      False
    )delim"
  );

  molecule.def(
    "assign_stereopermutator_randomly",
    [](Molecule& mol, AtomIndex a) {
      mol.assignStereopermutatorRandomly(a, randomnessEngine());
    },
    pybind11::arg("atom"),
    R"delim(
      Assigns an :class:`AtomStereopermutator` at random (assignments are
      weighted by relative statistical occurence).

      :param atom: Atom index of the stereopermutator to assign randomly.

      .. note::
         This function advances ``molassembler``'s global PRNG state.

      >>> # Assign an unspecified chiral center
      >>> mol = io.experimental.from_smiles("S[As](F)(Cl)(Br)(N)[H]")
      >>> as_index = 1
      >>> mol.stereopermutators.option(as_index).assigned is None
      True
      >>> mol.assign_stereopermutator_randomly(1)
      >>> mol.stereopermutators.option(as_index).assigned is None
      False
    )delim"
  );

  molecule.def(
    "assign_stereopermutator_randomly",
    [](Molecule& mol, const BondIndex& b) {
      mol.assignStereopermutatorRandomly(b, randomnessEngine());
    },
    pybind11::arg("bond_index"),
    R"delim(
      Assigns a :class:`BondStereopermutator` at random.

      :param bond_index: :class:`BondIndex` of the stereopermutator to assign randomly.

      .. note::
         This function advances ``molassembler``'s global PRNG state.

      >>> # Assign an unspecified double bond randomly
      >>> mol = io.experimental.from_smiles("CC=CC")
      >>> double_bond_index = BondIndex(1, 2)
      >>> assert mol.graph.bond_type(double_bond_index) == BondType.Double
      >>> mol.stereopermutators.option(double_bond_index).assigned is None
      True
      >>> mol.assign_stereopermutator_randomly(double_bond_index)
      >>> mol.stereopermutators.option(double_bond_index).assigned is None
      False
    )delim"
  );

  molecule.def(
    "canonicalize",
    &Molecule::canonicalize,
    pybind11::arg("components_bitmask") = AtomEnvironmentComponents::All,
    R"delim(
      Transform the molecule to a canonical form.

      :warning: Invalidates all external atom and bond indices.

      Molecule instances can be canonicalized. Graph canonicalization is an
      algorithm that reduces all isomorphic forms of an input graph into a
      canonical form. After canonicalization, isomorphism tests are reduced to
      mere identity tests.

      The canonicalization itself, however, is computationally at least as
      expensive as an isomorphism itself. Therefore, no expense is saved if an
      isomorphism test is to be computed only once for two molecules by
      canonizing both. Only if a molecule instance is to be a repeated
      candidate for isomorphism tests is there value in canonizing it.

      This library takes the approach of adding a tag to molecules that
      identifies which components of the graph and stereocenters have been used
      in the generation of the canonical form. This tag is voided with the use
      of any non-const member function. Pay close attention to the
      documentation of comparison member functions and operators to ensure that
      you are making good use of the provided shortcuts.

      Note that canonicalization information is only retained across IO
      boundaries using the JSON serialization variations.

      :param components_bitmask: The components of the molecular graph to
        include in the canonicalization procedure.
      :return: Flat index mapping/permutation from old indices to new

      >>> # Create two different representations of the same molecule
      >>> a = io.experimental.from_smiles("N[C@](Br)(O)C")
      >>> b = io.experimental.from_smiles("Br[C@](O)(N)C")
      >>> # a and be represent the same molecule, but have different vertex order
      >>> a == b # Equality operators perform an isomorphism for non-canonical pairs
      True
      >>> amap = a.canonicalize()
      >>> bmap = b.canonicalize()
      >>> amap == bmap # This shows the vertex order was different
      False
      >>> a == b # Equality operators perform a same-graph test for canonical pairs (faster)
      True
    )delim"
  );

  molecule.def(
    "remove_atom",
    &Molecule::removeAtom,
    pybind11::arg("atom"),
    R"delim(
      Remove an atom from the graph, including bonds to it, after checking
      that removing it is safe, i.e. the removal does not disconnect the graph.

      :warning: Invalidates all external atom and bond indices.

      :param atom: Atom to remove

      >>> m = Molecule() # Make H2
      >>> [a for a in m.graph.atoms()]
      [0, 1]
      >>> m.graph.can_remove(0) # We can remove a hydrogen from H2
      True
      >>> m.remove_atom(0)
      >>> m.graph.V  # We are left with just a hydrogen atom
      1
      >>> m.graph.E
      0
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

      :warning: Invalidates all external atom and bond indices.

      :param first_atom: First atom of the bond to be removed
      :param second_atom: Second atom of the bond to be removed

      >>> cyclopropane = io.experimental.from_smiles("C1CC1")
      >>> # In cyclopropane, we can remove a C-C bond without disconnecting the graph
      >>> cyclopropane.graph.can_remove(BondIndex(0, 1))
      True
      >>> V_before = cyclopropane.graph.V
      >>> E_before = cyclopropane.graph.E
      >>> cyclopropane.remove_bond(BondIndex(0, 1))
      >>> V_before - cyclopropane.graph.V # The number of atoms is unchanged
      0
      >>> E_before - cyclopropane.graph.E # We really only removed a bond
      1
      >>> # Note that now the valence of the carbon atoms where we removed
      >>> # the bond is... funky
      >>> cyclopropane.graph.degree(0)
      3
      >>> expected_bonds = [BondType.Single, BondType.Single, BondType.Single]
      >>> g = cyclopropane.graph
      >>> [g.bond_type(b) for b in g.bonds(0)] == expected_bonds
      True
      >>> cyclopropane.stereopermutators.option(0).shape == shapes.Shape.VacantTetrahedron
      True
    )delim"
  );

  molecule.def(
    "remove_bond",
    pybind11::overload_cast<const BondIndex&>(&Molecule::removeBond),
    pybind11::arg("bond_index"),
    R"delim(
      Remove a bond from the graph, after checking that removing it is safe,
      i.e. the removal does not disconnect the graph.

      :param bond_index: :class:`BondIndex` of the bond to be removed
    )delim"
  );

  molecule.def(
    "remove_permutator",
    &Molecule::removePermutator,
    pybind11::arg("bond_index"),
    R"delim(
      Remove a bond stereopermutator from the molecule, if present

      :param bond_index: Bond from which to remove the stereopermutator
      :returns: Whether a stereopermutator was removed
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
      :param bond_type: The new :class:`BondType`
      :return: Whether the bond already existed

      >>> # You really do have full freedom when it comes to your graphs:
      >>> h2 = Molecule()
      >>> _ = h2.set_bond_type(0, 1, BondType.Double) # Double bonded hydrogen atoms!
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

      >>> # Transform H2 into HF
      >>> import scine_utilities as utils
      >>> from copy import copy
      >>> H2 = Molecule()
      >>> HF = copy(H2)
      >>> HF.set_element_type(0, utils.ElementType.F)
      >>> HF == H2
      False
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

      >>> # Make methane square planar
      >>> from copy import copy
      >>> methane = io.experimental.from_smiles("C")
      >>> square_planar_methane = copy(methane)
      >>> square_planar_methane.set_shape_at_atom(0, shapes.Shape.Square)
      >>> methane == square_planar_methane
      False
    )delim"
  );

  molecule.def(
    "thermalize_stereopermutator",
    &Molecule::thermalizeStereopermutator,
    pybind11::arg("atom_index"),
    pybind11::arg("thermalization") = true,
    R"delim(
      Change the thermalization at an atom stereopermutator

      Alters the thermalization of stereopermutations at an atom
      stereopermutator.

      :param atom_index: Atom whose atom stereopermutator's thermalization to
        change
      :param thermalization: New status of thermalization to set
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
    R"delim(
      Read only access to the graph representation

      :rtype: :class:`Graph`
    )delim"
  );

  molecule.def_property_readonly(
    "stereopermutators",
    &Molecule::stereopermutators,
    R"delim(
      Read only access to the list of stereopermutators

      :rtype: :class:`StereopermutatorList`
    )delim"
  );

  molecule.def_property_readonly(
    "canonical_components",
    &Molecule::canonicalComponents,
    R"delim(
      Yields the components of the molecule that were used in a previous
      canonicalization. Can be ``None`` if the molecule was never
      canonicalized.

      :rtype: :class:`AtomEnvironmentComponents` or ``None``

      >>> # Canonicalize something and retrieve its canonical components
      >>> mol = io.experimental.from_smiles("C12(CCC1)COCC2")
      >>> mol.canonical_components is None
      True
      >>> _ = mol.canonicalize()
      >>> mol.canonical_components == AtomEnvironmentComponents.All
      True
    )delim"
  );

  molecule.def(
    "canonical_compare",
    &Molecule::canonicalCompare,
    pybind11::arg("other"),
    pybind11::arg("components_bitmask") = AtomEnvironmentComponents::All,
    R"delim(
      Modular comparison of this Molecule with another, assuming that both are
      in some (possibly partial) canonical form.

      For comparisons of fully canonical molecule pairs, regular equality
      comparison will just call this function with all environment components
      considered instead of performing a full isomorphism.

      This function is similar to modular_isomorphism, but faster, since if both
      molecules are in a canonical form, comparison does not require an
      isomorphism, but merely a same-graph test over the components used.

      :param other: The other (canonical) molecule to compare against
      :param components_bitmask: The components of an atom's environment to
        include in the comparison. You should use the same bitmask as when
        canonicalizing the molecules you are comparing here. It may be possible
        to use a bitmask with fewer components, but certainly not one with more.

      >>> # Bring two molecules into a partial canonical form and compare them
      >>> a = io.experimental.from_smiles("OCC")
      >>> b = io.experimental.from_smiles("SCC")
      >>> a == b
      False
      >>> # A and B are identical when considered purely by their graph
      >>> part = AtomEnvironmentComponents.Connectivity
      >>> _ = a.canonicalize(part)
      >>> _ = b.canonicalize(part)
      >>> a.canonical_compare(b, part)
      True
      >>> a == b # Partial canonicalization does not change the meaning of strict equality
      False
      >>> # Another pair that is identical save for a stereopermutation
      >>> c = io.experimental.from_smiles("N[C@](Br)(O)C")
      >>> d = io.experimental.from_smiles("N[C@@](Br)(O)C")
      >>> c == d # Strict equality includes stereopermutation
      False
      >>> part = AtomEnvironmentComponents.ElementsBondsAndShapes
      >>> _ = c.canonicalize(part)
      >>> _ = d.canonicalize(part)
      >>> c.canonical_compare(d, part) # Limited comparison yields equality
      True
    )delim"
  );

  molecule.def(
    "modular_isomorphism",
    &Molecule::modularIsomorphism,
    pybind11::arg("other"),
    pybind11::arg("components_bitmask"),
    R"delim(
      Modular comparison of this Molecule with another.

      This permits detailed specification of which elements of the molecular
      information you want to use in the comparison.

      Equality comparison is performed in several stages: First, at each atom
      position, a hash is computed that encompasses all local information that
      is specified to be used by the ``components_bitmask`` parameter. This
      hash is then used during graph isomorphism calculation to avoid finding
      an isomorphism that does not consider the specified factors.

      If an isomorphism is found, it is then validated. Bond orders and
      stereopermutators across both molecules are compared using the found
      isomorphism as an index map.

      Shortcuts to ``canonical_compare`` if ``components_bitmask`` matches the
      canonical components of both molecules (see ``canonical_components``).

      :param other: The molecule to compare against
      :param components_bitmask: The components of the molecule to use in the
        comparison

      :returns: None if the molecules are not isomorphic, a List[int] index
        mapping from self to other if the molecules are isomorphic.

      >>> a = io.experimental.from_smiles("OCC")
      >>> b = io.experimental.from_smiles("SCC")
      >>> a == b
      False
      >>> # A and B are identical when considered purely by their graph
      >>> a.modular_isomorphism(b, AtomEnvironmentComponents.Connectivity) is not None
      True
      >>> # Another pair that is identical save for a stereopermutation
      >>> c = io.experimental.from_smiles("N[C@](Br)(O)C")
      >>> d = io.experimental.from_smiles("N[C@@](Br)(O)C")
      >>> c == d # Strict equality includes stereopermutation
      False
      >>> c.modular_isomorphism(d, AtomEnvironmentComponents.ElementsBondsAndShapes) is not None
      True
    )delim"
  );

  molecule.def(pybind11::self == pybind11::self);
  molecule.def(pybind11::self != pybind11::self);

  /* Integration with IPython / Jupyter */
  if(graphvizInPath()) {
    molecule.def(
      "_repr_svg_",
      &::pipeSVG,
      "Generates an SVG representation of the molecule"
    );
  }

  // Shell integration
  molecule.def(
    "__str__",
    [](const Molecule& mol) -> std::string {
      std::string repr = "Molecule of elemental composition ";
      repr += Scine::Utils::generateChemicalFormula(mol.graph().elementCollection());
      const unsigned A = mol.stereopermutators().A();
      const unsigned B = mol.stereopermutators().B();

      if(A > 0) {
        repr += " with " + std::to_string(A) + " atom ";

        if(B > 0) {
          repr += "and " + std::to_string(B) + " bond ";
        }
        repr += "stereopermutator";
        if(A + B > 1) {
          repr += "s";
        }
      }

      return repr;
    },
    "Generate a string summary of the molecule components"
  );

  /* Direct copying support */
  molecule.def(
    "__copy__",
    [](const Molecule& mol) -> Molecule {
      return Molecule(mol);
    }
  );
  molecule.def(
    "__deepcopy__",
    [](const Molecule& mol, pybind11::dict /* memo */) -> Molecule {
      return Molecule(mol);
    },
    pybind11::arg("memo")
  );

  /* Pickling support */
  molecule.def(
    pybind11::pickle(
      [](const Molecule& mol) {
        JsonSerialization serialization(mol);
        std::string serialized = serialization;
        return serialized;
      },
      [](const std::string& serialized) {
        JsonSerialization serialization(serialized);
        Molecule mol = serialization;
        return mol;
      }
    )
  );
}
