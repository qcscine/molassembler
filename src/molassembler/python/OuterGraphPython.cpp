/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "molassembler/Cycles.h"
#include "molassembler/Graph.h"

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/FormulaGenerator.h"
#include "Utils/Typenames.h"

void init_outer_graph(pybind11::module& m) {
  using namespace Scine::molassembler;
  pybind11::class_<Graph> outerGraph(
    m,
    "Graph",
    R"delim(
      Molecular graph in which atoms are vertices and bonds are edges.

      >>> import scine_utilities as utils
      >>> ethane = io.experimental.from_smiles("CC")
      >>> g = ethane.graph
      >>> g.atoms_of_element(utils.ElementType.C)
      [0, 1]
      >>> g.degree(0)
      4
      >>> g.can_remove(0)
      False
      >>> g.can_remove(BondIndex(0, 1))
      False
      >>> hydrogen_indices = g.atoms_of_element(utils.ElementType.H)
      >>> can_remove = lambda a : g.can_remove(a)
      >>> all(map(can_remove, hydrogen_indices))
      True
      >>> g.N
      8
      >>> g.B
      7
    )delim"
  );

  outerGraph.def(
    "adjacent",
    &Graph::adjacent,
    pybind11::arg("first_atom"),
    pybind11::arg("second_atom"),
    R"delim(
      Returns whether two atoms are bonded

      >>> ethane = io.experimental.from_smiles("CC")
      >>> ethane.graph.degree(0)
      4
      >>> [ethane.graph.adjacent(0, a) for a in range(1, ethane.graph.N)]
      [True, True, True, True, False, False, False]
    )delim"
  );

  outerGraph.def(
    "atoms_of_element",
    &Graph::atomsOfElement,
    pybind11::arg("element_type"),
    R"delim(
      Returns atoms matching an element type

      >>> import scine_utilities as utils
      >>> ethanol = io.experimental.from_smiles("CCO")
      >>> ethanol.graph.atoms_of_element(utils.ElementType.O)
      [2]
      >>> ethanol.graph.atoms_of_element(utils.ElementType.C)
      [0, 1]
    )delim"
  );

  outerGraph.def(
    "bond_orders",
    &Graph::bondOrders,
    R"delim(
      Generates a BondOrderCollection representation of the molecule connectivity

      >>> # Convert acetaldehyde's graph into a floating point bond order matrix
      >>> import scine_utilities as utils
      >>> acetaldehyde = io.experimental.from_smiles("CC=O")
      >>> bo = acetaldehyde.graph.bond_orders()
      >>> bo.empty()
      False
      >>> bo.get_order(0, 1) # The order between the carbon atoms
      1.0
      >>> bo.get_order(1, 2) # The order between a carbon and oxygen
      2.0
    )delim"
  );

  outerGraph.def(
    "bond_type",
    &Graph::bondType,
    pybind11::arg("bond_index"),
    R"delim(
      Fetches the :class:`BondType` at a particular :class:`BondIndex`

      >>> # Look at some bond orders of an interesting model compound
      >>> compound = io.experimental.from_smiles("[Co]1(C#N)(C#O)C=C1")
      >>> compound.graph.bond_type(BondIndex(0, 1)) # Co-CN bond
      BondType.Single
      >>> compound.graph.bond_type(BondIndex(0, 5)) # Co-C=C bond
      BondType.Eta
      >>> compound.graph.bond_type(BondIndex(5, 6)) # C=C bond
      BondType.Double
      >>> compound.graph[BondIndex(1, 2)] # C#N bond by bond subsetting
      BondType.Triple
    )delim"
  );

  outerGraph.def(
    "can_remove",
    pybind11::overload_cast<AtomIndex>(&Graph::canRemove, pybind11::const_),
    pybind11::arg("atom"),
    R"delim(
      Returns whether an atom can be removed without disconnecting the graph

      >>> # In graph terms, articulation vertices cannot be removed
      >>> methane = io.experimental.from_smiles("C")
      >>> methane.graph.can_remove(0) # We cannot remove the central carbon
      False
      >>> all([methane.graph.can_remove(i) for i in range(1, 5)]) # But hydrogens!
      True
    )delim"
  );

  outerGraph.def(
    "can_remove",
    pybind11::overload_cast<const BondIndex&>(&Graph::canRemove, pybind11::const_),
    pybind11::arg("bond_index"),
    R"delim(
      Returns whether a bond can be removed without disconnecting the graph

      >>> # In graph terms, bridge edges cannot be removed
      >>> import scine_utilities as utils
      >>> from itertools import combinations
      >>> cyclopropane = io.experimental.from_smiles("C1CC1")
      >>> carbon_atoms = cyclopropane.graph.atoms_of_element(utils.ElementType.C)
      >>> cc_bonds = [BondIndex(a, b) for (a, b) in combinations(carbon_atoms, 2)]
      >>> can_remove = lambda b: cyclopropane.graph.can_remove(b)
      >>> all(map(can_remove, cc_bonds)) # We can remove any one of the bonds
      True
      >>> cyclopropane.remove_bond(cc_bonds[0]) # Remove one C-C bond
      >>> any(map(can_remove, cc_bonds[1:])) # Can we still remove any of the others?
      False
    )delim"
  );

  outerGraph.def_property_readonly(
    "cycles",
    &Graph::cycles,
    "Fetch a reference to cycles information"
  );

  outerGraph.def(
    "degree",
    &Graph::degree,
    pybind11::arg("atom"),
    R"delim(
      Returns the number of bonds incident upon an atom.

      >>> # A silly example
      >>> model = io.experimental.from_smiles("CNO[H]")
      >>> [model.graph.degree(i) for i in range(0, 4)]
      [4, 3, 2, 1]
    )delim"
  );

  // TODO rename to elements
  outerGraph.def(
    "element_collection",
    &Graph::elementCollection,
    R"delim(
      Generates an ElementCollection representation of the molecule's atoms' element types

      >>> # Some isotopes
      >>> import scine_utilities as utils
      >>> m = io.experimental.from_smiles("[1H]C([2H])([3H])[H]")
      >>> m.graph.element_collection()
      [ElementType.H1, ElementType.C, ElementType.D, ElementType.T, ElementType.H]
    )delim"
  );

  outerGraph.def(
    "element_type",
    &Graph::elementType,
    pybind11::arg("atom"),
    R"delim(
      Fetch the element type of an atom

      >>> # Some isotopes
      >>> import scine_utilities as utils
      >>> m = io.experimental.from_smiles("[1H]C([2H])([3H])[H]")
      >>> m.graph.element_type(0)
      ElementType.H1
      >>> m.graph.element_type(2)
      ElementType.D
      >>> m.graph[4] # Subsettable wih atom indices to get element types
      ElementType.H
    )delim"
  );

  outerGraph.def_property_readonly("N", &Graph::N, "The number of atoms in the graph");
  outerGraph.def_property_readonly("B", &Graph::B, "The number of bonds in the graph");

  outerGraph.def(
    "split_along_bridge",
    &Graph::splitAlongBridge,
    pybind11::arg("bridge_bond"),
    R"delim(
      Determine which atoms belong to either side of a bond

      >>> # Hypothetically splitting a model compound
      >>> m = io.experimental.from_smiles("CN")
      >>> m.graph.split_along_bridge(BondIndex(0, 1))
      ([0, 2, 3, 4], [1, 5, 6])
    )delim"
  );

  outerGraph.def(
    "atoms",
    [](const Graph& graph) {
      auto atomsRange = graph.atoms();
      return pybind11::make_iterator(
        std::move(atomsRange.first),
        std::move(atomsRange.second)
      );
    },
    R"delim(
      Iterate through all valid atom indices of the graph

      Fully equivalent to: ``range(graph.N)``
    )delim"
  );

  outerGraph.def(
    "bonds",
    [](const Graph& graph) {
      auto bondsRange = graph.bonds();
      return pybind11::make_iterator(
        std::move(bondsRange.first),
        std::move(bondsRange.second)
      );
    },
    R"delim(
      Iterate through all valid bond indices of the graph

      >>> import scine_utilities as utils
      >>> model = io.experimental.from_smiles("F/C=C/I")
      >>> [b for b in model.graph.bonds()]
      [(0, 1), (1, 2), (2, 3), (1, 4), (2, 5)]
    )delim"
  );

  outerGraph.def(
    "adjacents",
    [](const Graph& graph, AtomIndex a) {
      auto adjacentsRange = graph.adjacents(a);
      return pybind11::make_iterator(
        std::move(adjacentsRange.first),
        std::move(adjacentsRange.second)
      );
    },
    pybind11::arg("a"),
    R"delim(
      Iterate through all adjacent atom indices of an atom

      >>> import scine_utilities as utils
      >>> m = io.experimental.from_smiles("NC")
      >>> [a for a in m.graph.adjacents(0)]
      [1, 2, 3]
      >>> element = lambda a: m.graph.element_type(a)
      >>> [element(a) for a in m.graph.adjacents(0)]
      [ElementType.C, ElementType.H, ElementType.H]
    )delim"
  );

  outerGraph.def(
    "bonds",
    [](const Graph& graph, AtomIndex a) {
      auto bondsRange = graph.bonds(a);
      return pybind11::make_iterator(
        std::move(bondsRange.first),
        std::move(bondsRange.second)
      );
    },
    pybind11::arg("a"),
    R"delim(
      Iterate through all incident bonds of an atom

      >>> import scine_utilities as utils
      >>> m = io.experimental.from_smiles("NC")
      >>> [b for b in m.graph.bonds(0)]
      [(0, 1), (0, 2), (0, 3)]
    )delim"
  );

  outerGraph.def(
    "__repr__",
    [](const Graph& graph) {
      return (
        "Graph of elemental composition "
        + Scine::Utils::generateChemicalFormula(graph.elementCollection())
      );
    }
  );

  outerGraph.def(
    "__getitem__",
    [](const Graph& g, const AtomIndex i) -> Scine::Utils::ElementType {
      return g.elementType(i);
    }
  );

  outerGraph.def(
    "__getitem__",
    [](const Graph& g, const BondIndex i) -> BondType {
      return g.bondType(i);
    }
  );
}
