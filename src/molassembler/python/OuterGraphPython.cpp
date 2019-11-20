/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "molassembler/Cycles.h"
#include "molassembler/OuterGraph.h"

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Geometry/FormulaGenerator.h"
#include "Utils/Typenames.h"

void init_outer_graph(pybind11::module& m) {
  using namespace Scine::molassembler;
  pybind11::class_<OuterGraph> outerGraph(
    m,
    "Graph",
    R"delim(
      Molecular graph in which atoms are vertices and bonds are edges.

      >>> import molassembler as masm
      >>> import scine_utils_os as utils
      >>> ethane = masm.io.experimental.from_smiles("CC")
      >>> g = ethane.graph
      >>> g.atoms_of_element(utils.ElementType.C)
      [0, 1]
      >>> g.degree(0)
      4
      >>> g.can_remove(0)
      False
      >>> g.can_remove(masm.BondIndex(0, 1))
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
    &OuterGraph::adjacent,
    pybind11::arg("first_atom"),
    pybind11::arg("second_atom"),
    R"delim(
      Returns whether two atoms are bonded

      >>> import molassembler as masm
      >>> ethane = masm.io.experimental.from_smiles("CC")
      >>> ethane.graph.degree(0)
      4
      >>> [ethane.graph.adjacent(0, a) for a in range(1, ethane.graph.N)]
      [True, True, True, True, False, False, False]
    )delim"
  );

  outerGraph.def(
    "atoms_of_element",
    &OuterGraph::atomsOfElement,
    pybind11::arg("element_type"),
    R"delim(
      Returns atoms matching an element type

      >>> import molassembler as masm
      >>> import scine_utils_os as utils
      >>> ethanol = masm.io.experimental.from_smiles("CCO")
      >>> ethanol.graph.atoms_of_element(utils.ElementType.O)
      [2]
      >>> ethanol.graph.atoms_of_element(utils.ElementType.C)
      [0, 1]
    )delim"
  );

  outerGraph.def(
    "bond_orders",
    &OuterGraph::bondOrders,
    R"delim(
      Generates a BondOrderCollection representation of the molecule connectivity

      >>> # Convert acetaldehyde's graph into a floating point bond order matrix
      >>> import molassembler as masm
      >>> import scine_utils_os as utils
      >>> acetaldehyde = masm.io.experimental.from_smiles("CC=O")
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
    &OuterGraph::bondType,
    pybind11::arg("bond_index"),
    R"delim(
      Fetches the :class:`BondType` at a particular :class:`BondIndex`

      >>> # Look at some bond orders of an interesting model compound
      >>> import molassembler as masm
      >>> compound = masm.io.experimental.from_smiles("[Co]1(C#N)(C#O)C=C1")
      >>> compound.graph.bond_type(masm.BondIndex(0, 1)) # Co-CN bond
      BondType.Single
      >>> compound.graph.bond_type(masm.BondIndex(0, 5)) # Co-C=C bond
      BondType.Eta
      >>> compound.graph.bond_type(masm.BondIndex(5, 6)) # C=C bond
      BondType.Double
      >>> compound.graph.bond_type(masm.BondIndex(1, 2)) # C#N bond
      BondType.Triple
    )delim"
  );

  outerGraph.def(
    "can_remove",
    pybind11::overload_cast<AtomIndex>(&OuterGraph::canRemove, pybind11::const_),
    pybind11::arg("atom"),
    "Returns whether an atom can be removed without disconnecting the graph"
  );

  outerGraph.def(
    "can_remove",
    pybind11::overload_cast<const BondIndex&>(&OuterGraph::canRemove, pybind11::const_),
    pybind11::arg("bond_index"),
    "Returns whether a bond can be removed without disconnecting the graph"
  );

  outerGraph.def_property_readonly(
    "cycles",
    &OuterGraph::cycles,
    "Fetch a reference to cycles information"
  );

  outerGraph.def(
    "degree",
    &OuterGraph::degree,
    pybind11::arg("atom"),
    "Returns the number of bonds incident upon an atom"
  );

  outerGraph.def(
    "element_collection",
    &OuterGraph::elementCollection,
    "Generates an ElementCollection representation of the molecule's atoms' element types"
  );

  outerGraph.def(
    "element_type",
    &OuterGraph::elementType,
    pybind11::arg("atom"),
    "Fetch the element type of an atom"
  );

  outerGraph.def_property_readonly("N", &OuterGraph::N, "The number of atoms in the graph");
  outerGraph.def_property_readonly("B", &OuterGraph::B, "The number of bonds in the graph");

  outerGraph.def(
    "split_along_bridge",
    &OuterGraph::splitAlongBridge,
    pybind11::arg("bridge_bond"),
    "Determine which atoms belong to either side of a bond"
  );

  outerGraph.def(
    "atoms",
    [](const OuterGraph& graph) {
      auto atomsRange = graph.atoms();
      return pybind11::make_iterator(
        std::move(atomsRange.first),
        std::move(atomsRange.second)
      );
    },
    "Iterate through all valid atom indices of the graph"
  );

  outerGraph.def(
    "bonds",
    [](const OuterGraph& graph) {
      auto bondsRange = graph.bonds();
      return pybind11::make_iterator(
        std::move(bondsRange.first),
        std::move(bondsRange.second)
      );
    },
    "Iterate through all valid bond indices of the graph"
  );

  outerGraph.def(
    "adjacents",
    [](const OuterGraph& graph, AtomIndex a) {
      auto adjacentsRange = graph.adjacents(a);
      return pybind11::make_iterator(
        std::move(adjacentsRange.first),
        std::move(adjacentsRange.second)
      );
    },
    pybind11::arg("a"),
    "Iterate through all adjacent atom indices of an atom"
  );

  outerGraph.def(
    "bonds",
    [](const OuterGraph& graph, AtomIndex a) {
      auto bondsRange = graph.bonds(a);
      return pybind11::make_iterator(
        std::move(bondsRange.first),
        std::move(bondsRange.second)
      );
    },
    pybind11::arg("a"),
    "Iterate through all incident bonds of an atom"
  );

  outerGraph.def(
    "__repr__",
    [](const OuterGraph& graph) {
      return (
        "Graph of elemental composition "
        + Scine::Utils::generateChemicalFormula(graph.elementCollection())
      );
    }
  );
}
