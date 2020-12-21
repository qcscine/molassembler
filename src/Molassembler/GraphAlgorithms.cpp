/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/GraphAlgorithms.h"

#include "Molassembler/Graph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Graph/EditDistance.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Subgraphs.h"
#include "Molassembler/Molecule/MolGraphWriter.h"
#include "Molassembler/Graph/McSplit.h"

#include "Utils/Geometry/ElementInfo.h"

#include "boost/graph/graphviz.hpp"
#include "boost/process/child.hpp"
#include "boost/process/io.hpp"
#include "boost/process/search_path.hpp"

/* GraphAlgorithms.h is the public API point for algorithms that act on
 * Graph. In Graph/GraphAlgorithms.h, those algorithms are implemented on
 * the PrivateGraph, which is not accessible using the public API.
 */

namespace Scine {
namespace Molassembler {
namespace {

using UnorderedIndexMap = std::unordered_map<AtomIndex, MinimalReactionEdits::ComponentIndexPair>;
std::pair<PrivateGraph, UnorderedIndexMap> condense(const GraphList& list) {
  std::pair<PrivateGraph, UnorderedIndexMap> p;
  unsigned component = 0;
  for(const auto& graphRef : list) {
    const auto indexMap = p.first.merge(graphRef.get().inner());
    for(auto pair : indexMap) {
      p.second.emplace(pair.second, std::make_pair(component, pair.first));
    }
    ++component;
  }
  assert(p.second.size() == p.first.V());
  return p;
}

MinimalGraphEdits minimalEdits(
  const PrivateGraph& a,
  const PrivateGraph& b,
  const EditCost& cost,
  const bool preconditioningSubgraphConnected
) {
  /* Precondition the graph edit distance algorithm with maximum common
   * subgraph matches
   *
   * TODO The McSplit doesn't seem to offer all mcs (just one), and in
   * principle we need to try all of them and minimize the distance over them.
   */
  constexpr bool mcsLabelEdges = true;
  const auto preconditioning = GraphAlgorithms::McSplit::mcs(
    a,
    b,
    preconditioningSubgraphConnected,
    mcsLabelEdges
  );
  Subgraphs::IndexMap precondition;
  for(auto p : preconditioning) {
    precondition.insert(Subgraphs::IndexMap::value_type(p.first, p.second));
  }

  // Template for minimization over subgraphs:
  // const auto edits = Temple::accumulate(
  //   Subgraphs::maximum(lhs.first, rhs.first),
  //   Edits { std::numeric_limits<unsigned>::max(), {}},
  //   [&](Edits carry, const auto& preconditioning) -> MinimalGraphEdits {
  //     std::cout << "MCS precondition size: " << preconditioning.size() << "\n";
  //     EditForest forest {a, b, cost, precondition};
  //     const unsigned distance = forest.g[forest.result].costSum;

  //     if(distance < carry.distance) {
  //       return preconditionedEdits;
  //     }

  //     return carry;
  //   }
  // );

  using EditForest = GraphAlgorithms::EditDistanceForest;
  EditForest forest {a, b, cost, precondition};

  static_assert(
    EditForest::epsilon == MinimalGraphEdits::epsilon,
    "Epsilon values must match between EditDistanceForest and MinimalGraphEdits"
  );
  constexpr auto epsilon =  EditForest::epsilon;

  MinimalGraphEdits edits;
  edits.distance = forest.g[forest.result].costSum;
  const auto traversal = forest.traverse(forest.result);
  edits.indexMap = traversal.bVertices;
  std::reverse(std::begin(edits.indexMap), std::end(edits.indexMap));

  // Compute the full list of non-zero cost edits
  const auto sizes = std::minmax({a.V(), b.V()});
  const unsigned edgeAlterationCost = cost.edgeAlteration();
  // Vertices from a
  for(AtomIndex i = 0; i < sizes.first; ++i) {
    const AtomIndex j = edits.indexMap[i];
    if(j == epsilon && cost.vertexAlteration() > 0) {
      edits.vertexEdits.emplace_back(i, j, cost.vertexAlteration());
    } else {
      const unsigned elementCost = cost.elementSubstitution(a.elementType(i), b.elementType(j));
      if(elementCost > 0) {
        edits.vertexEdits.emplace_back(i, j, elementCost);
      }
    }

    // Implied edges
    for(AtomIndex k = 0; k < i; ++k) {
      const AtomIndex l = edits.indexMap.at(k);
      const auto aBond = EditForest::bondTypeOption(i, k, a);
      const auto bBond = EditForest::bondTypeOption(j, l, b);
      if(aBond && bBond) {
        const unsigned edgeCost = cost.bondSubstitution(*aBond, *bBond);
        if(edgeCost > 0) {
          edits.edgeEdits.emplace_back(BondIndex {i, k}, BondIndex {j, l}, edgeCost);
        }
      } else if(static_cast<bool>(aBond) ^ static_cast<bool>(bBond)) {
        if(edgeAlterationCost > 0) {
          edits.edgeEdits.emplace_back(BondIndex {i, k}, BondIndex {j, l}, edgeAlterationCost);
        }
      }
    }
  }

  // Extra vertices from b (epsilon -> j)
  if(edgeAlterationCost != 0) {
    for(AtomIndex i = sizes.first; i < sizes.second; ++i) {
      const AtomIndex j = edits.indexMap[i];
      for(AtomIndex k = 0; k < i; ++k) {
        const AtomIndex l = edits.indexMap[k];
        if(l == epsilon) {
          continue;
        }
        if(auto bEdgeOption = b.edgeOption(j, l)) {
          edits.edgeEdits.emplace_back(
            BondIndex {epsilon, k < sizes.first ? k : epsilon},
            BondIndex {j, l},
            edgeAlterationCost
          );
        }
      }
    }
  }

  return edits;
}

struct ReactionGraphWriter : public MolGraphWriter {
  using ComponentBondIndex = MinimalReactionEdits::EdgeEdit::ComponentBondIndex;

  ReactionGraphWriter(
    const PrivateGraph* graph,
    UnorderedIndexMap passIdxMap,
    const std::vector<ComponentBondIndex>& edits
  ) : MolGraphWriter(graph, nullptr),
      idxMap(std::move(passIdxMap)),
      edgeEdits(std::begin(edits), std::end(edits)) {}

  std::string vertexLabel(const PrivateGraph::Vertex v) const override {
    const Utils::ElementType e = graphPtr->elementType(v);
    const auto w = idxMap.at(v).second;

    // Do not mark element type for hydrogen or carbon, those are commonplace
    if(e == Utils::ElementType::H || e == Utils::ElementType::C) {
      return std::to_string(w);
    }

    const std::string symbolString = Utils::ElementInfo::symbol(e);
    return symbolString + std::to_string(w);
  }

  std::string edgeColor(const PrivateGraph::Edge& e) const override {
    const BondIndex b {graphPtr->source(e), graphPtr->target(e)};

    const ComponentBondIndex maybeEdited {
      idxMap.at(b.first),
      idxMap.at(b.second)
    };

    if(edgeEdits.count(maybeEdited) > 0) {
      return "tomato";
    }

    return "black";
  }

  UnorderedIndexMap idxMap;
  std::set<ComponentBondIndex> edgeEdits;
};

template<typename T, typename U>
void extractBoostStream(T& t, U& u) {
#if BOOST_VERSION >= 107000
  /* NOTE: This implementation of buffer transfers in boost process has a bug
   * that isn't fixed before Boost 1.70.
   */
  t << u.rdbuf();
#else
  // Workaround: cast to a parent class implementing rdbuf() correctly.
  using BasicIOSReference = std::basic_ios<char, std::char_traits<char>>&;
  // Feed the results into our ostream
  t << static_cast<BasicIOSReference>(u).rdbuf();
#endif
}

} // namespace

std::vector<unsigned> distance(AtomIndex source, const Graph& graph) {
  if(source > graph.V()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return GraphAlgorithms::distance(source, graph.inner());
}

PredecessorMap shortestPaths(AtomIndex source, const Graph& graph) {
  if(source > graph.V()) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  return PredecessorMap {
    GraphAlgorithms::shortestPaths(source, graph.inner())
  };
}

std::vector<AtomIndex> PredecessorMap::path(const AtomIndex target) const {
  std::vector<AtomIndex> pathVertices;
  AtomIndex position = target;
  while(predecessors.at(position) != position) {
    pathVertices.push_back(position);
    position = predecessors.at(position);
  }
  pathVertices.push_back(position);

  std::reverse(
    std::begin(pathVertices),
    std::end(pathVertices)
  );

  return pathVertices;
}

constexpr AtomIndex MinimalGraphEdits::epsilon;

MinimalGraphEdits minimalEdits(const Graph& a, const Graph& b, const EditCost& cost) {
  constexpr bool preconditioningSubgraphConnected = true;
  return minimalEdits(a.inner(), b.inner(), cost, preconditioningSubgraphConnected);
}

MinimalReactionEdits reactionEdits(const GraphList& lhsGraphs, const GraphList& rhsGraphs) {
  const auto lhs = condense(lhsGraphs);
  const auto rhs = condense(rhsGraphs);

  // Check preconditions
  auto lhsElements = lhs.first.elementCollection();
  auto rhsElements = rhs.first.elementCollection();
  if(lhsElements.size() != rhsElements.size()) {
    throw std::logic_error("Unequal atom count in reaction. Playing at alchemy?");
  }
  Temple::sort(lhsElements);
  Temple::sort(rhsElements);
  if(lhsElements != rhsElements) {
    throw std::logic_error("Element composition of reaction sides different. Playing at alchemy?");
  }

  const auto edits = minimalEdits(
    lhs.first,
    rhs.first,
    ElementsConservedCost {},
    false
  );

  MinimalReactionEdits multi;
  multi.distance = edits.distance;

  const unsigned V = edits.indexMap.size();
  assert(V == lhs.first.V());
  for(unsigned i = 0; i < V; ++i) {
    // NOTE: This just assumes there are no vertex epsilons due to the high cost
    multi.indexMap.emplace(
      lhs.second.at(i),
      rhs.second.at(edits.indexMap.at(i))
    );
  }

  multi.vertexEdits = Temple::map(
    edits.vertexEdits,
    [&](const auto& edit) -> MinimalReactionEdits::VertexEdit {
      return {
        lhs.second.at(edit.first),
        rhs.second.at(edit.second),
        edit.cost
      };
    }
  );

  multi.edgeEdits = Temple::map(
    edits.edgeEdits,
    [&](const auto& edit) -> MinimalReactionEdits::EdgeEdit {
      return {
        {lhs.second.at(edit.first.first), lhs.second.at(edit.first.second)},
        {rhs.second.at(edit.second.first), rhs.second.at(edit.second.second)},
        edit.cost
      };
    }
  );

  return multi;
}

std::string reactionGraphvizSvg(
  const GraphList& lhsGraphs,
  const GraphList& rhsGraphs,
  const MinimalReactionEdits& edits
) {
  namespace bp = boost::process;
  const auto dot = bp::search_path("dot");
  if(dot.empty()) {
    throw std::runtime_error("dot graphviz binary not found in PATH");
  }

  // Preconditions: Find the graphviz binaries we need
  const auto gvpack = bp::search_path("gvpack");
  if(gvpack.empty()) {
    throw std::runtime_error("gvpack graphviz binary not found in PATH");
  }

  const auto neato = bp::search_path("neato");
  if(neato.empty()) {
    throw std::runtime_error("neato graphviz binary not found in PATH");
  }

  const auto lhs = condense(lhsGraphs);
  const auto rhs = condense(rhsGraphs);

  const auto lhsEdgeEdits = Temple::map(edits.edgeEdits,
    [](const auto& edit) { return edit.first; }
  );
  ReactionGraphWriter lhsWriter(&lhs.first, lhs.second, lhsEdgeEdits);

  const auto rhsEdgeEdits = Temple::map(edits.edgeEdits,
    [](const auto& edit) { return edit.second; }
  );
  ReactionGraphWriter rhsWriter(&rhs.first, rhs.second, rhsEdgeEdits);

  /* Now for the piping. Need to imitate the following sequence:
   *
   *   cat <(dot lhs.dot) <(dot rhs.dot)
   *   | gvpack -array_uc1 -Glayout=neato 2>/dev/null
   *   | neato -n2 -Tsvg
   *
   */
  bp::opstream graphs;
  auto dotGraph = [&](auto& g, auto& writer) {
    bp::opstream graphInput;
    bp::ipstream dotOutput;
    boost::write_graphviz(
      graphInput,
      g,
      writer,
      writer,
      writer
    );

    bp::child dotChild(
      dot,
      bp::std_in < graphInput,
      bp::std_out > dotOutput,
      bp::std_err > bp::null
    );

    graphInput.flush();
    graphInput.pipe().close();
    dotChild.wait();

    // Reeeally inefficient, but I don't know what I'm doing
    std::stringstream tmp;
    extractBoostStream(tmp, dotOutput);
    graphs << tmp.str();
  };

  dotGraph(lhs.first.bgl(), lhsWriter);
  dotGraph(rhs.first.bgl(), rhsWriter);

  bp::ipstream packed_o;
  bp::opstream packed_i;
  bp::ipstream svg;

  bp::child pack(
    gvpack.string() + " -array_uc1 -Glayout=neato",
    bp::std_in < graphs,
    bp::std_out > packed_o,
    bp::std_err > bp::null
  );
  graphs.flush();
  graphs.pipe().close();
  pack.wait();

  extractBoostStream(packed_i, packed_o);
  bp::child render(
    neato.string() + " -n2 -Tsvg",
    bp::std_in < packed_i,
    bp::std_out > svg,
    bp::std_err > bp::null
  );
  packed_i.flush();
  packed_i.pipe().close();
  render.wait();

  std::stringstream os;
  extractBoostStream(os, svg);
  return os.str();
}

} // namespace Molassembler
} // namespace Scine
