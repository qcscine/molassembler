/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/IO/SmilesEmitter.h"

#include "Molassembler/IO/SmilesCommon.h"
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

AtomIndex spanningTreeRoot(const Molecule& mol) {
  using IndexDistancePair = std::pair<AtomIndex, unsigned>;
  return Temple::accumulate(
    mol.graph().atoms(),
    IndexDistancePair {0, 0},
    [&](const IndexDistancePair& carry, const AtomIndex i) -> IndexDistancePair {
      const auto distances = distance(i, mol.graph());
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

struct Emitter {
  struct RingClosure {
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

  static AtomIndex root(const Molecule& mol) {
    using IndexDistancePair = std::pair<AtomIndex, unsigned>;
    return Temple::accumulate(
      mol.graph().atoms(),
      IndexDistancePair {0, 0},
      [&](const IndexDistancePair& carry, const AtomIndex i) -> IndexDistancePair {
        const auto distances = distance(i, mol.graph());
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

  static bool organicAliphatic(
    const Utils::ElementType e,
    const unsigned hydrogenCount
  ) {
    if(!isValenceFillElement(e)) {
      return false;
    }
  }

  Emitter(const Molecule& mol) : molecule(mol) {
    const PrivateGraph& inner = mol.graph().inner();

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

    AtomIndex suggestedRoot = root(mol);
    const auto findIter = subsumedHydrogenAtoms.find(suggestedRoot);
    if(findIter != std::end(subsumedHydrogenAtoms)) {
      suggestedRoot = findIter->second;
    }

    const auto distances = distance(suggestedRoot, mol.graph());
    const AtomIndex partner = std::max_element(
      std::begin(distances),
      std::end(distances)
    ) - std::begin(distances);
    longestPath = shortestPaths(suggestedRoot, mol.graph()).path(partner);

    const PrivateGraph::Vertex V = inner.V();
    std::vector<PrivateGraph::Vertex> predecessors (V);
    boost::prim_minimum_spanning_tree(
      inner.bgl(),
      boost::make_iterator_property_map(
        predecessors.begin(),
        boost::get(boost::vertex_index, inner.bgl())
      ),
      boost::root_vertex(suggestedRoot)
    );

    // Can't DFS a predecessor map - need to make an actual tree instead
    spanning = SpanningTree(V);
    // Add edges from predecessor list
    for(PrivateGraph::Vertex v = 0; v < V; ++v) {
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
  }

  void addClosure(const AtomIndex u, const AtomIndex v) {
    auto listIter = closures.find(u);
    if(listIter == closures.end()) {
      listIter = closures.emplace_hint(listIter, u);
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

    if(!isIsotope) {
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

  void dfs(const AtomIndex v, std::vector<AtomIndex>::iterator& pathIter) {
    // Write the atom symbol for v
    smiles += atomSymbol(v);
    /* Are we on the longest path? If so, impose special ordering on dfs
     * that descends along the longest path last.
     */
    if(v == *pathIter) {
      ++pathIter;
      const AtomIndex nextOnLongestPath = *pathIter;
      auto iters = boost::out_edges(v, spanning);
      for(; iters.first != iters.second; ++iters.first) {
        const AtomIndex w = boost::target(*iters.first, spanning);
        smiles += "(";
        if(w != nextOnLongestPath) {
          // TODO Write edge type
          dfs(w, pathIter);
        }
        smiles += ")";
      }
      // TODO write edge type
      dfs(nextOnLongestPath, pathIter);
    } else {
      // TODO need to know which is last so we can add parentheses to all priors
      auto iters = boost::out_edges(v, spanning);
      for(; iters.first != iters.second; ++iters.first) {
        const AtomIndex w = boost::target(*iters.first, spanning);
        // TODO Write edge type
        dfs(w, pathIter);
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

  std::string smiles;
  SpanningTree spanning;
  std::vector<AtomIndex> longestPath;
  std::unordered_set<AtomIndex> aromaticAtoms;
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
