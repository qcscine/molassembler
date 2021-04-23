/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Utils/Pybind.h"

#include "Molassembler/Graph.h"
#include "Molassembler/GraphAlgorithms.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Temple/Functional.h"

#include <fstream>

using namespace Scine;
using namespace Molassembler;

void init_graph_algorithms(pybind11::module& m) {
  m.def(
    "distance",
    &distance,
    pybind11::arg("source"),
    pybind11::arg("graph"),
    R"delim(
      Calculates graph distances from a single atom index to all others

      >>> m = io.experimental.from_smiles("CC(CC)C")
      >>> distances = distance(1, m.graph)
    )delim"
  );

  m.def(
    "sites",
    [](const Graph& graph, const AtomIndex v) {
      return GraphAlgorithms::sites(graph.inner(), v);
    },
    pybind11::arg("graph"),
    pybind11::arg("atom"),
    R"delim(
      Returns adjacents of an atom of the graph grouped into sites

      Sites consisting of multiple atoms are haptic.
    )delim"
  );

  pybind11::class_<PredecessorMap> predecessor_map(m, "PredecessorMap");

  predecessor_map.def_readonly("predecessors", &PredecessorMap::predecessors);
  predecessor_map.def(
    "path",
    &PredecessorMap::path,
    pybind11::arg("target"),
    R"delim(
      Generate path vertices to target vertex

      :param target: Target vertex to generate path to
      :returns: Vertex path starting at source and including target
    )delim"
  );

  m.def(
    "shortest_paths",
    &shortestPaths,
    pybind11::arg("source"),
    pybind11::arg("graph"),
    R"delim(
      Generate predecessor map containing shortest paths to each vertex in a graph

      :param source: Source vertex to generate shortest paths from
      :param graph: Graph in which to generate shortest paths

      >>> m = io.experimental.from_smiles("CC(CC)C")
      >>> predecessors = shortest_paths(1, m.graph)
      >>> predecessors.path(0)
      [1, 0]
      >>> predecessors.path(3)
      [1, 2, 3]
    )delim"
  );

  // TODO documentation here, and in C++ parts
  pybind11::class_<EditCost> editCost(
    m,
    "EditCost",
    R"delim(
      Cost functor for graph edit distance calculations.
    )delim"
  );

  pybind11::class_<FuzzyCost, EditCost> fuzzyCost(
    m,
    "FuzzyCost",
    R"delim(
      Cost functor for fuzzy graph edit distances. All costs (vertex and edge
      alteration, element substitution and bond order substitution) are
      unitary.
    )delim"
  );
  fuzzyCost.def(pybind11::init<>());

  pybind11::class_<ElementsConservedCost, EditCost> elementsConservedCost(
    m,
    "ElementsConservedCost",
    R"delim(
      Cost functor for graph edit distances conserving element types. Edge
      alteration and bond order substitution costs are unitary, but vertex
      alteration and element substitution are of cost 100.
    )delim"
  );
  elementsConservedCost.def(pybind11::init<>());

  pybind11::class_<MinimalGraphEdits> editsCls(
    m,
    "MinimalGraphEdits",
    R"delim(
      Data class containing the results of a minimal graph edit distance
      calculation.
    )delim"
  );
  editsCls.def_readonly(
    "distance",
    &MinimalGraphEdits::distance,
    R"delim(
      Total cost of the minimal edits. Symmetric under ordering changes of the
      graph input to the calculation.
    )delim"
  );
  editsCls.def_readonly(
    "index_map",
    &MinimalGraphEdits::indexMap,
    R"delim(
      Flat index mapping from the left graph to the right graph. Can contain
      epsilon values indicating vertex deletion. Values at indices starting at
      the size of the first graph are inserted vertices.
    )delim"
  );
  editsCls.def_property_readonly(
    "vertex_edits",
    [](const MinimalGraphEdits& edits) {
      return Temple::map(
        edits.vertexEdits,
        [](const auto& edit) {
          return std::make_tuple(edit.first, edit.second, edit.cost);
        }
      );
    },
    R"delim(
      Non-zero cost vertex edits. Tuple of the left vertex index, the right
      vertex index, and the incurred cost. May contain epsilon values.
    )delim"
  );
  editsCls.def_property_readonly(
    "edge_edits",
    [](const MinimalGraphEdits& edits) {
      return Temple::map(
        edits.edgeEdits,
        [](const auto& edit) {
          return std::make_tuple(edit.first, edit.second, edit.cost);
        }
      );
    },
    R"delim(
      Non-zero cost edge edits. Tuple of the left bond index, the right bond
      index, and the incurred cost. May contain epsilon values.
    )delim"
  );
  editsCls.def_readonly_static(
    "epsilon",
    &MinimalGraphEdits::epsilon,
    R"delim(
      Sentinel value indicating vertex deletion in the index map or in bond
      indices of the edits.
    )delim"
  );
  m.def(
    "minimal_edits",
    &minimalEdits,
    pybind11::arg("a"),
    pybind11::arg("b"),
    Utils::Arg("cost") = FuzzyCost {},
    R"delim(
      Minimal graph edits

      Calculates a minimal set of vertex or edge insertions/deletions and
      their incurred cost to edit one graph into another.

      The applied cost function encourages fuzzy matching while preferring
      label substitution over insertions and deletions.

      A maximum common connected subgraph is used to precondition the exact graph
      edit distance algorithm. Otherwise, the exact algorithm quickly becomes
      intractable past 10 vertices due to combinatorial space explosion and
      rapid memory exhaustion.

      :param a: First graph to calculate edit distance for
      :param b: Second graph to calculate edit distance for
      :param cost: Cost function for minimal edits. Defaults to FuzzyCost.
        Alternative is ElementsConservedCost.
    )delim"
  );

  pybind11::class_<MinimalReactionEdits> reactionEditsCls(
    m,
    "MinimalReactionEdits",
    R"delim(
      Data class for multiple-graph minimal edits in reactions
    )delim"
  );
  reactionEditsCls.def_readonly(
    "distance",
    &MinimalReactionEdits::distance,
    R"delim(
      Total cost of the minimal edits. Symmetric under ordering changes of the
      graph input to the calculation.
    )delim"
  );
  reactionEditsCls.def_readonly(
    "index_map",
    &MinimalReactionEdits::indexMap,
    R"delim(
      Index map from the left graph to the right graph. Each side of the mapping
      is composed of a component index (indicating which of the graphs of the
      side of the reaction the mapped vertex is part of) and its atom index.
    )delim"
  );
  reactionEditsCls.def_property_readonly(
    "vertex_edits",
    [](const MinimalReactionEdits& reactionEdits) {
      return Temple::map(
        reactionEdits.vertexEdits,
        [](const auto& edit) {
          return std::make_tuple(edit.first, edit.second, edit.cost);
        }
      );
    },
    R"delim(
      Non-zero cost vertex edits. Tuple of left and right component and vertex
      index pairs, and the incurred cost.
    )delim"
  );
  reactionEditsCls.def_property_readonly(
    "edge_edits",
    [](const MinimalReactionEdits& reactionEdits) {
      return Temple::map(
        reactionEdits.edgeEdits,
        [](const auto& edit) {
          return std::make_tuple(edit.first, edit.second, edit.cost);
        }
      );
    },
    R"delim(
      Non-zero cost edge edits. Tuple of left and right component and vertex
      bond index pairs, and the incurred cost.
    )delim"
  );

  m.def(
    "reaction_edits",
    &reactionEdits,
    R"delim(
      Minimal reaction edits

      Calculates a minimal set of vertex or edge insertions/deletions and
      their incurred cost to edit two sides of a reaction into one another.

      The applied cost function heavily discourages element type label
      substitutions and vertex insertions or deletions.

      A maximum common subgraph is used to precondition the exact graph
      edit distance algorithm. Otherwise, the exact algorithm quickly becomes
      intractable past 10 vertices due to combinatorial space explosion and
      rapid memory exhaustion.

      :param lhs: List of graphs of the left side of the reaction
      :param rhs: List of graphs of the right side of the reaction

      :raises RuntimeError: If the number of atoms in both sides is unequal or
        the element composition of both sides is different.

      :returns: distance, index mapping and non-zero cost edit lists
    )delim"
  );

  struct ReactionEditSvg { std::string svg; };
  pybind11::class_<ReactionEditSvg> svgCls(m, "ReactionEditSvg");
  svgCls.def_property_readonly("svg", [](const ReactionEditSvg& c) { return c.svg; });
  svgCls.def("_repr_svg_", [](const ReactionEditSvg& c) { return c.svg; });
  svgCls.def(
    "write",
    [](const ReactionEditSvg& c, const std::string& fname) {
      std::ofstream file(fname);
      file << c.svg;
    }
  );

  m.def(
    "reaction_edits_svg",
    [](
      const GraphList& lhs,
      const GraphList& rhs,
      const MinimalReactionEdits& edits
    ) -> ReactionEditSvg {
      return { reactionGraphvizSvg(lhs, rhs, edits) };
    },
    pybind11::arg("lhs"),
    pybind11::arg("rhs"),
    pybind11::arg("reaction_edits"),
    R"delim(
      Generate a graphviz representation of changes in a chemical reaction

      Requires the graphviz binaries dot, neato and gvpack to be available in
      the PATH.

      :raises RuntimeError: If the required graphviz binaries are not found.
    )delim"
  );
}
