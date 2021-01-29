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
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/Shapes/Properties.h"

#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Permutations.h"
#include "Utils/Geometry/ElementInfo.h"
#include "boost/graph/prim_minimum_spanning_tree.hpp"

// #include <fstream>
// #include "boost/graph/graphviz.hpp"
// #include "Molassembler/Temple/Stringify.h"
// #include <iostream>

#include <unordered_set>

namespace Scine {
namespace Molassembler {
namespace IO {
namespace Experimental {
namespace {

struct Unity {
  using G = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::bidirectionalS
  >;
  using key_type = typename G::edge_descriptor;
  using value_type = double;
  using reference = double;
  using category = boost::readable_property_map_tag;
};

inline double get(const Unity& /* u */, const Unity::key_type& /* e */) {
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
  const unsigned inverseOrder = 7 - order;

  /* Weight edges by the graph degree (of heavy atoms only) to incentivize ring
   * closures at fulcrums of the heavy atom graph. Important for normalization
   * of e.g. biphenyl.
   */
  auto heavyDegree = [&](const AtomIndex i) -> unsigned {
    return Temple::accumulate(
      g.adjacents(i),
      0U,
      [&](const unsigned carry, const AtomIndex j) -> unsigned {
        if(
          g.elementType(j) == Utils::ElementType::H
          && g.degree(j) == 1
        ) {
          return carry;
        }

        return carry + 1;
      }
    );
  };

  return inverseOrder + heavyDegree(g.source(e)) + heavyDegree(g.target(e));
}

const std::unordered_map<unsigned, unsigned>& smilesStereodescriptorShapeMap(
  const Shapes::Shape shape
) {
  switch(shape) {
    case Shapes::Shape::Tetrahedron: {
      static std::unordered_map<unsigned, unsigned> tetrahedronMap {
        {0, 1},
        {1, 2},
      };
      return tetrahedronMap;
    }
    case Shapes::Shape::Square: {
      static std::unordered_map<unsigned, unsigned> squareMap {
        {0, 1},
        {3, 2},
        {22, 3}
      };
      return squareMap;
    }
    case Shapes::Shape::TrigonalBipyramid: {
      static std::unordered_map<unsigned, unsigned> trigBipyramidMap {
        {32, 1},
        {38, 2},
        {34, 3},
        {44, 4},
        {40, 5},
        {46, 6},
        {64, 7},
        {70, 8},
        {8, 9},
        {10, 10},
        {14, 11},
        {20, 12},
        {16, 13},
        {22, 14},
        {2, 15},
        {4, 16},
        {0, 17},
        {6, 18},
        {18, 19},
        {12, 20},
      };
      return trigBipyramidMap;
    }
    case Shapes::Shape::Octahedron: {
      static std::unordered_map<unsigned, unsigned> octahedronMap {
        {152, 1},
        {566, 2},
        {154, 3},
        {158, 4},
        {164, 5},
        {160, 6},
        {166, 7},
        {176, 8},
        {178, 9},
        {542, 10},
        {662, 11},
        {202, 12},
        {668, 13},
        {446, 14},
        {470, 15},
        {686, 16},
        {590, 17},
        {710, 18},
        {184, 19},
        {190, 20},
        {208, 21},
        {692, 22},
        {596, 23},
        {716, 24},
        {304, 25},
        {310, 26},
        {328, 27},
        {694, 28},
        {598, 29},
        {718, 30},
      };
      return octahedronMap;
    }
    default: throw std::out_of_range("No shape map for selected shape");
  }
}


unsigned smilesStereodescriptor(
  const Shapes::Shape shape,
  const std::vector<Shapes::Vertex>& order
) {
  const auto allRotations = Shapes::Properties::generateAllRotations(
    shape,
    order
  );
  std::unordered_set<unsigned> permutationIndices;
  for(const auto& indices : allRotations) {
    permutationIndices.insert(Temple::permutationIndex(indices));
  }

  for(const auto& pair : smilesStereodescriptorShapeMap(shape)) {
    if(permutationIndices.count(pair.first) > 0) {
      return pair.second;
    }
  }

  throw std::out_of_range("No smiles stereodescriptor found!");
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

  static bool heteroatom(const Utils::ElementType e) {
    return e != Utils::ElementType::H && e != Utils::ElementType::C;
  }

  Emitter(const Molecule& mol) : molecule(mol) {
    const PrivateGraph& inner = mol.graph().inner();

    /* Find terminal hydrogen atoms that can be rewritten into heavy atom
     * labels, implicit or explicit
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
    AtomIndex suggestedRoot = root(mol);
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

    AtomIndex partner = std::max_element(
      std::begin(distances),
      std::end(distances)
    ) - std::begin(distances);
    longestPath = PredecessorMap {predecessors}.path(partner);

    // Smiles normalization: Begin on a heteroatom, if possible
    if(
      heteroatom(molecule.graph().elementType(partner))
      && !heteroatom(molecule.graph().elementType(suggestedRoot))
    ) {
      std::swap(suggestedRoot, partner);
      std::reverse(std::begin(longestPath), std::end(longestPath));
    }

    /* Now for the second minimum spanning tree calculation to figure out the
     * directionality of the tree
     */
    boost::prim_minimum_spanning_tree(
      spanning,
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

    markAromaticAtoms();
  }

  void markAromaticAtoms() {
    const std::unordered_set<Shapes::Shape> flatShapes {
      Shapes::Shape::Bent,
      Shapes::Shape::EquilateralTriangle
    };

    for(const auto& cycleEdges : molecule.graph().cycles()) {
      const bool permutatorsEverywhere = Temple::all_of(
        cycleEdges,
        [&](const BondIndex& bond) -> bool {
          return static_cast<bool>(molecule.stereopermutators().option(bond));
        }
      );
      if(!permutatorsEverywhere) {
        continue;
      }
      auto cycleAtoms = makeRingIndexSequence(cycleEdges);
      const bool flat = Temple::all_of(
        cycleAtoms,
        [&](const AtomIndex i) -> bool {
          auto permutator = molecule.stereopermutators().option(i);
          if(!permutator) {
            return false;
          }
          return flatShapes.count(permutator->getShape()) > 0;
        }
      );
      if(!flat) {
        continue;
      }

      for(const AtomIndex i : cycleAtoms) {
        aromaticAtoms.insert(i);
      }
      for(const BondIndex& b : cycleEdges) {
        aromaticBonds.insert(b);
      }
    }
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

  bool isSmilesStereogenic(const AtomIndex v) const {
    auto permutator = molecule.stereopermutators().option(v);
    if(!permutator) {
      return false;
    }
    switch(permutator->getShape()) {
      // These are the only shapes that smiles has descriptors for
      case Shapes::Shape::Tetrahedron:
      case Shapes::Shape::Square:
      case Shapes::Shape::TrigonalBipyramid:
      case Shapes::Shape::Octahedron:
        return permutator->numAssignments() > 1;
      default:
        return false;
    }
  }

  std::string atomSymbol(
    const AtomIndex v,
    const std::vector<AtomIndex>& descendantOrder
  ) const {
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

    const bool isStereogenic = isSmilesStereogenic(v);

    if(!isIsotope && !isStereogenic && isValenceFillElement(e)) {
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

          const BondType type = molecule.graph().bondType(BondIndex {v, i});
          const int order = Bond::bondOrderMap.at(static_cast<unsigned>(type));
          return count + order;
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

    // Bracket atom
    std::string elementStr = "[";
    if(isIsotope) {
      elementStr += std::to_string(Utils::ElementInfo::A(e));
    }
    elementStr += Utils::ElementInfo::symbol(e);
    if(isStereogenic) {
      elementStr += atomStereodescriptor(v, descendantOrder);
    }
    if(hydrogenCount > 1) {
      elementStr += "H" + std::to_string(hydrogenCount);
    } else if(hydrogenCount == 1) {
      elementStr += "H";
    }
    elementStr += "]";
    return elementStr;
  }

  std::string atomStereodescriptor(
    const AtomIndex v,
    const std::vector<AtomIndex>& descendantOrder
  ) const {
    unsigned hydrogenCount = 0;
    const auto countIter = subsumedHydrogenCount.find(v);
    if(countIter != subsumedHydrogenCount.end()) {
      hydrogenCount = countIter->second;
    }

    auto permutator = molecule.stereopermutators().option(v);
    std::vector<SiteIndex> siteSequence;
    /* - Use the shape map to place the substituent order index of the smiles
     *   at the shape vertices (via the site indices), taking extra care of
     *   hydrogen counts (which are ordered first in the smiles sequence)
     * - Generate all rotations of that permutation and find a smiles
     *   stereodescriptor for that shape that matches
     */
    if(hydrogenCount > 0) {
      SiteIndex i {0};
      for(const auto& siteIndices : permutator->getRanking().sites) {
        if(siteIndices.size() == 1) {
          const AtomIndex siteAtom = siteIndices.front();
          if(
            molecule.graph().elementType(siteAtom) == Utils::ElementType::H
            && molecule.graph().degree(siteAtom) == 1
          ) {
            siteSequence.push_back(i);
          }
        }
        ++i;
      }
    }
    for(auto iters = boost::in_edges(v, spanning); iters.first != iters.second; ++iters.first) {
      const AtomIndex parent = boost::source(*iters.first, spanning);
      siteSequence.push_back(permutator->getRanking().getSiteIndexOf(parent));
    }
    for(const AtomIndex i : descendantOrder) {
      siteSequence.push_back(permutator->getRanking().getSiteIndexOf(i));
    }

    const auto orderAtShapeVertices = Temple::map(
      siteSequence,
      [&](const SiteIndex i) -> Shapes::Vertex {
        return permutator->getShapePositionMap().at(i);
      }
    );

    const unsigned stereodescriptor = smilesStereodescriptor(
      permutator->getShape(),
      Shapes::Properties::inverseRotation(orderAtShapeVertices)
    );

    switch(permutator->getShape()) {
      case Shapes::Shape::Tetrahedron: {
        if(stereodescriptor == 1) {
          return "@";
        }
        return "@@";
      }
      case Shapes::Shape::Square: {
        return "@SP" + std::to_string(stereodescriptor);
      }
      case Shapes::Shape::TrigonalBipyramid: {
        return "@TB" + std::to_string(stereodescriptor);
      }
      case Shapes::Shape::Octahedron: {
        return "@OH" + std::to_string(stereodescriptor);
      }
      default: std::exit(1); // Unreachable by precondition of isStereogenic
    }
  }

  void writeEdgeType(BondIndex b) {
    switch(molecule.bondType(b)) {
      case BondType::Single: {
        /* Single bond order isn't generally written, except if it's a bond
         * between aromatic atoms that isn't aromatic itself, e.g. in biphenyl
         * c1ccccc1-c2ccccc2
         */
        if(
          aromaticAtoms.count(b.first) > 0
          && aromaticAtoms.count(b.second) > 0
          && aromaticBonds.count(b) == 0
        ) {
          smiles += "-";
        }
        break;
      }
      case BondType::Double: { smiles += "="; break; }
      case BondType::Triple: { smiles += "#"; break; }
      // NOTE: Eta and bond orders four through six aren't explicitly written
      default: break;
    }
  }

  void dfs(const AtomIndex v, std::vector<AtomIndex>::iterator& pathIter) {
    // Determine the order of the substituents indicated by the MST
    std::vector<AtomIndex> descendantOrder;
    for(auto iters = boost::out_edges(v, spanning); iters.first != iters.second; ++iters.first) {
      descendantOrder.push_back(boost::target(*iters.first, spanning));
    }
    const bool onMainBranch = (v == *pathIter);
    if(onMainBranch) {
      // Impose requirement that the main branch is descended onto last
      ++pathIter;
      if(pathIter != std::end(longestPath)) {
        const auto nextMainBranchAtom = std::find(
          std::rbegin(descendantOrder),
          std::rend(descendantOrder),
          *pathIter
        );
        std::iter_swap(nextMainBranchAtom, std::rbegin(descendantOrder));
      }
    }

    // Write the atom symbol for v
    smiles += atomSymbol(v, descendantOrder);
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

    // Descend into substituents
    for(const AtomIndex w : descendantOrder) {
      const bool last = (w == descendantOrder.back());
      if(!last) {
        smiles += "(";
      }
      writeEdgeType(BondIndex {v, w});
      dfs(w, pathIter);
      if(!last) {
        smiles += ")";
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
  std::unordered_set<BondIndex, boost::hash<BondIndex>> aromaticBonds;
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
