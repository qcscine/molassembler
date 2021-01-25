/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/IO/SmilesEmitter.h"

#include "Molassembler/IO/SmilesCommon.h"
#include "Molassembler/Modeling/BondDistance.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/GraphAlgorithms.h"

#include "Molassembler/Temple/Functional.h"
#include "Utils/Geometry/ElementInfo.h"
#include "boost/graph/prim_minimum_spanning_tree.hpp"

#include <unordered_set>

namespace Scine {
namespace Molassembler {
namespace IO {
namespace Experimental {
namespace {

struct Unity {
  using key_type = PrivateGraph::Edge;
  using value_type = double;
  using reference = double;
  using category = boost::readable_property_map_tag;
};

inline double get(const Unity& /* u */, const PrivateGraph::Edge& /* e */) {
  return 1.0;
}

struct InverseBondOrder {
  using key_type = PrivateGraph::Edge;
  using value_type = double;
  using reference = double;
  using category = boost::readable_property_map_tag;

  explicit InverseBondOrder(const PrivateGraph& g) : graphRef(g) {}

  std::reference_wrapper<const PrivateGraph> graphRef;
};

inline double get(const InverseBondOrder& u, const PrivateGraph::Edge& e) {
  const PrivateGraph& g = u.graphRef.get();
  const BondType type = g.bondType(e);
  const double order = Bond::bondOrderMap.at(static_cast<unsigned>(type));
  return 7 - order;
}

} // namespace

struct Emitter {
  struct RingClosure {
    RingClosure() = default;
    explicit RingClosure(AtomIndex i) : partner(i) {}

    // Bond partner of ring-closing bond
    AtomIndex partner;
    // Marked number in depth-first traversal
    boost::optional<unsigned> number;
  };

  using ClosureMap = std::unordered_map<AtomIndex, std::vector<RingClosure>>;

  using SpanningTree = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::bidirectionalS
  >;

  Emitter(const Molecule& mol) : molecule(mol) {
    const PrivateGraph& inner = mol.graph().inner();

    /* TODO
     * - The procedure here is faulty! Need to break ring-closing bonds first,
     *   then figure out the longest path in the resulting graph.
     * - This might be achieved by
     *   - MST on the original graph,
     *   - construct a copy of the graph with all edges in the MST in
     *     bidirectional form
     *   - determine the longest path in the graph
     *   - MST with the root of that longest path
     *   - Reduce edges to the directionality of the MST
     *
     * - Some SMILES normalization issues that could be addressed
     *   - Avoid making ring-closures on double or triple bonds
     *     - Maybe weight edges during MST generation by inverse order (single
     *       is 6, double is 5 etc.). That should encourage marking low-order
     *       bonds as ring-closures over higher-order bonds.
     *   - Start on a heteroatom if possible (OCC over CCO)
     *     - Reverse the longest path if necessary
     */

    /* Figure out the rewriting of hydrogen atoms into heavy atoms labels,
     * implicit or explicit
     */
    std::unordered_map<AtomIndex, AtomIndex> subsumedHydrogenAtoms;
    for(const AtomIndex i : inner.vertices()) {
      if(
        inner.elementType(i) == Utils::ElementType::H
        && inner.degree(i) == 1
      ) {
        for(const AtomIndex j : inner.adjacents(i)) {
          // Cannot subsume the terminal hydrogen in another hydrogen atom
          if(inner.elementType(j) != Utils::ElementType::H) {
            subsumedHydrogenAtoms.emplace(i, j);
          }
        }
      }
    }

    // Rewrite subsumption mapping into hydrogen counts at the target vertex
    for(const auto& iterPair : subsumedHydrogenAtoms) {
      auto findIter = subsumedHydrogenCount.find(iterPair.second);
      if(findIter == subsumedHydrogenCount.end()) {
        subsumedHydrogenCount.emplace_hint(findIter, iterPair.second, 1);
      } else {
        ++findIter->second;
      }
    }

    /* TODO maybe need to place the root vertex of this MST at the 'centroid'
     * of the graph, i.e. the vertex with the lowest average distance to all
     * other vertices? Unsure of the influence of this arbitrary choice here.
     */
    const PrivateGraph::Vertex V = inner.V();
    std::vector<PrivateGraph::Vertex> predecessors (V);
    boost::prim_minimum_spanning_tree(
      inner.bgl(),
      boost::make_iterator_property_map(
        predecessors.begin(),
        boost::get(boost::vertex_index, inner.bgl())
      ),
      boost::root_vertex(0).
      weight_map(InverseBondOrder {inner})
    );

    spanning = SpanningTree(V);

    // Add edges from predecessor list, but bidirectional!
    for(const PrivateGraph::Vertex v : inner.vertices()) {
      const PrivateGraph::Vertex predecessor = predecessors.at(v);

      // Skip unreachable vertices and the root
      if(predecessor == v) {
        continue;
      }

      // Skip subsumed hydrogen atoms
      if(subsumedHydrogenAtoms.count(v) > 0) {
        continue;
      }

      boost::add_edge(predecessor, v, spanning);
      boost::add_edge(v, predecessor, spanning);
    }

    /* Figure out ring-closing bonds from predecessor list: If neither of the
     * bonds vertices has the other atom as its predecessor, it is not in the
     * minimum spanning tree and represents a ring closure
     */
    for(const auto& bond : inner.edges()) {
      const PrivateGraph::Vertex u = inner.source(bond);
      const PrivateGraph::Vertex v = inner.target(bond);
      if(predecessors.at(u) != v && predecessors.at(v) != u) {
        addClosure(u, v);
        addClosure(v, u);
      }
    }

    // Figure out the longest path in the acyclic graph
    const AtomIndex suggestedRoot = root(mol);
    std::fill(std::begin(predecessors), std::end(predecessors), suggestedRoot);
    std::vector<unsigned> distances(inner.V(), 0);
    boost::breadth_first_search(
      spanning,
      suggestedRoot,
      boost::visitor(
        boost::make_bfs_visitor(
          std::make_pair(
            boost::record_distances(&distances[0], boost::on_tree_edge {}),
            boost::record_predecessors(&predecessors[0], boost::on_tree_edge {})
          )
        )
      )
    );

    const AtomIndex partner = std::max_element(
      std::begin(distances),
      std::end(distances)
    ) - std::begin(distances);
    longestPath = PredecessorMap {predecessors}.path(partner);

    /* Now for the second minimum spanning tree calculation to figure out the
     * directionality of the tree
     */
    boost::prim_minimum_spanning_tree(
      inner.bgl(),
      boost::make_iterator_property_map(
        predecessors.begin(),
        boost::get(boost::vertex_index, inner.bgl())
      ),
      boost::root_vertex(suggestedRoot).
      weight_map(Unity {})
    );

    // Impose directionality of the MST on the graph
    for(PrivateGraph::Vertex v = 0; v < V; ++v) {
      const PrivateGraph::Vertex predecessor = predecessors.at(v);
      // Skip unreachable vertices and the root
      if(predecessor == v) {
        continue;
      }

      boost::remove_edge(v, predecessor, spanning);
    }

    // TODO Figure out aromatic atoms and bonds
  }

  std::vector<unsigned> distance(AtomIndex a) const {
    std::vector<unsigned> distances(boost::num_vertices(spanning), 0);

    boost::breadth_first_search(
      spanning,
      a,
      boost::visitor(
        boost::make_bfs_visitor(
          boost::record_distances(&distances[0], boost::on_tree_edge {})
        )
      )
    );

    return distances;
  }

  AtomIndex root(const Molecule& mol) const {
    using IndexDistancePair = std::pair<AtomIndex, unsigned>;
    return Temple::accumulate(
      mol.graph().atoms(),
      IndexDistancePair {0, 0},
      [&](const IndexDistancePair& carry, const AtomIndex i) -> IndexDistancePair {
        const auto distances = distance(i);
        const unsigned maxDistance = *std::max_element(
          std::begin(distances),
          std::end(distances)
        );

        if(maxDistance > carry.second) {
          return std::make_pair(i, maxDistance);
        }

        return carry;
      }
    ).first;
  }

  void addClosure(const AtomIndex u, const AtomIndex v) {
    auto listIter = closures.find(u);
    if(listIter == closures.end()) {
      listIter = closures.emplace_hint(listIter, u, 0);
    }
    listIter->second.emplace_back(v);
  }

  std::string atomSymbol(const AtomIndex v) const {
    /* Decide whether to emit either
     * - an aliphatic organic symbol
     * - an aromatic organic symbol
     * - a bracket atom
     *   - Necessary if isotopic or non-standard hydrogen count
     */
    const Utils::ElementType e = molecule.elementType(v);
    const bool isIsotope = Utils::ElementInfo::base(e) != e;
    unsigned hydrogenCount = 0;
    const auto countIter = subsumedHydrogenCount.find(v);
    if(countIter != subsumedHydrogenCount.end()) {
      hydrogenCount = countIter->second;
    }

    if(!isIsotope && isValenceFillElement(e)) {
      const bool isAromatic = aromaticAtoms.count(v) > 0;

      const int valence = Temple::accumulate(
        molecule.graph().adjacents(v),
        0U,
        [&](const int count, const AtomIndex i) -> int {
          if(
            molecule.graph().elementType(i) == Utils::ElementType::H
            && molecule.graph().degree(i) == 1
          ) {
            return count;
          }

          return count + 1;
        }
      );

      const unsigned valenceFill = valenceFillElementImplicitHydrogenCount(
        valence,
        e,
        isAromatic
      );

      if(valenceFill == hydrogenCount) {
        std::string symbol = Utils::ElementInfo::symbol(e);
        if(isAromatic) {
          std::transform(
            std::begin(symbol),
            std::end(symbol),
            std::begin(symbol),
            [](unsigned char c) { return std::tolower(c); }
          );
        }
        return symbol;
      }
    }

    // Bracket atom fallback
    std::string elementStr = "[";
    if(isIsotope) {
      elementStr += std::to_string(Utils::ElementInfo::A(e));
    }
    elementStr += Utils::ElementInfo::symbol(e);
    if(hydrogenCount > 0) {
      elementStr += "H" + std::to_string(hydrogenCount);
    }
    elementStr += "]";
    return elementStr;
  }

  void writeEdgeType(BondIndex b) {
    /* TODO special case of single bonds between aromatic atoms if the bond
     * isn't aromatic, e.g. in biphenyl c1ccccc1-c2ccccc2
     */
    const BondType type = molecule.bondType(b);
    if(type == BondType::Double) {
      smiles += "=";
    }
    if(type == BondType::Triple) {
      smiles += "#";
    }
  }

  void dfs(const AtomIndex v, std::vector<AtomIndex>::iterator& pathIter) {
    // Write the atom symbol for v
    smiles += atomSymbol(v);
    // Ring closing bonds
    auto closuresIter = closures.find(v);
    if(closuresIter != closures.end()) {
      for(RingClosure& closure : closuresIter->second) {
        if(!closure.number) {
          closure.number = closureIndex;

          // Find matching closure for partner
          const auto partnerClosuresIter = closures.find(closure.partner);
          assert(partnerClosuresIter != closures.end());
          const auto partnerIter = std::find_if(
            std::begin(partnerClosuresIter->second),
            std::end(partnerClosuresIter->second),
            [&](const RingClosure& hay) {
              return hay.partner == v;
            }
          );
          assert(partnerIter != std::end(partnerClosuresIter->second));
          assert(!partnerIter->number);

          // Mark partner number
          partnerIter->number = closureIndex;

          ++closureIndex;
        }

        assert(closure.number);
        const unsigned number = closure.number.value();
        if(number < 10) {
          smiles += std::to_string(number);
        } else {
          smiles += "%" + std::to_string(number);
        }
      }
    }

    /* Are we on the longest path? If so, impose special ordering on dfs
     * that descends along the longest path last.
     */
    if(v == *pathIter) {
      ++pathIter;
      if(pathIter == std::end(longestPath)) {
        return;
      }
      const AtomIndex nextOnLongestPath = *pathIter;
      auto iters = boost::out_edges(v, spanning);
      for(; iters.first != iters.second; ++iters.first) {
        const AtomIndex w = boost::target(*iters.first, spanning);
        if(w != nextOnLongestPath) {
          smiles += "(";
          writeEdgeType(BondIndex {v, w});
          dfs(w, pathIter);
          smiles += ")";
        }
      }
      writeEdgeType(BondIndex {v, nextOnLongestPath});
      dfs(nextOnLongestPath, pathIter);
    } else {
      auto iters = boost::out_edges(v, spanning);
      for(; iters.first != iters.second; ++iters.first) {
        const bool last = iters.first + 1 != iters.second;
        if(!last) {
          smiles += "(";
        }

        const AtomIndex w = boost::target(*iters.first, spanning);
        writeEdgeType(BondIndex {v, w});
        dfs(w, pathIter);
        if(!last) {
          smiles += ")";
        }
      }
    }
  }

  std::string dfs() {
    /* Depth-first search with custom ordering that descends into the longest
     * path last at each vertex
     */
    auto pathIter = std::begin(longestPath);
    dfs(longestPath.front(), pathIter);
    return smiles;
  }

  unsigned closureIndex = 1;
  std::string smiles;
  SpanningTree spanning;
  std::vector<AtomIndex> longestPath;
  std::unordered_set<AtomIndex> aromaticAtoms;
  // TODO std::unordered_set<BondIndex> aromaticBonds;
  ClosureMap closures;
  std::unordered_map<AtomIndex, unsigned> subsumedHydrogenCount;
  const Molecule& molecule;
};

std::string emitSmiles(const Molecule& mol) {
  Emitter emitter(mol);
  return emitter.dfs();
}

} // namespace Experimental
} // namespace IO
} // namespace Molassembler
} // namespace Scine
