/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"

#include "molassembler/Cycles.h"
#include "molassembler/OuterGraph.h"

#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Typenames.h"

void init_outer_graph(pybind11::module& m) {
  using namespace Scine::molassembler;
  pybind11::class_<OuterGraph> outerGraph(m, "Graph", "Molecular graph");

  outerGraph.def("adjacent", &OuterGraph::adjacent, "Returns whether two atoms are bonded");
  outerGraph.def("bondOrders", &OuterGraph::bondOrders, "Generates a BondOrderCollection representation of the molecule connectivity");
  outerGraph.def("bondType", &OuterGraph::bondType, "Fetches the bond type at a particular bond");
  outerGraph.def(
    "canRemove",
    pybind11::overload_cast<AtomIndex>(&OuterGraph::canRemove, pybind11::const_),
    "Returns whether an atom can be removed without disconnecting the graph"
  );
  outerGraph.def(
    "canRemove",
    pybind11::overload_cast<const BondIndex&>(&OuterGraph::canRemove, pybind11::const_),
    "Returns whether a bond can be removed without disconnecting the graph"
  );
  outerGraph.def_property_readonly("cycles", &OuterGraph::cycles, "Fetch a reference to cycles information");
  outerGraph.def("degree", &OuterGraph::degree, "Returns the number of bonds incident upon an atom");
  outerGraph.def("elementCollection", &OuterGraph::elementCollection, "Generates an ElementCollection representation of the molecule's atoms' element types");
  outerGraph.def("elementType", &OuterGraph::elementType, "Fetch the element type of an atom");
  outerGraph.def("N", &OuterGraph::N, "The number of atoms in the graph");
  outerGraph.def("B", &OuterGraph::B, "The number of bonds in the graph");

  outerGraph.def(
    "atoms",
    [](const OuterGraph& graph) {
      auto atomsRange = graph.atoms();
      return pybind11::make_iterator(
        std::move(atomsRange.first),
        std::move(atomsRange.second)
      );
    }
  );

  outerGraph.def(
    "bonds",
    [](const OuterGraph& graph) {
      auto bondsRange = graph.bonds();
      return pybind11::make_iterator(
        std::move(bondsRange.first),
        std::move(bondsRange.second)
      );
    }
  );

  outerGraph.def(
    "adjacents",
    [](const OuterGraph& graph, AtomIndex a) {
      auto adjacentsRange = graph.adjacents(a);
      return pybind11::make_iterator(
        std::move(adjacentsRange.first),
        std::move(adjacentsRange.second)
      );
    }
  );

  outerGraph.def(
    "bonds",
    [](const OuterGraph& graph, AtomIndex a) {
      auto bondsRange = graph.bonds(a);
      return pybind11::make_iterator(
        std::move(bondsRange.first),
        std::move(bondsRange.second)
      );
    }
  );
}
