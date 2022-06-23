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

#include <fstream>
#include "boost/graph/graphviz.hpp"
#include "Molassembler/Temple/Stringify.h"
#include <iostream>

#include <unordered_set>

namespace Scine {
namespace Molassembler {
namespace IO {
namespace Experimental {

struct RingClosure {
  RingClosure() = default;
  explicit RingClosure(AtomIndex i) : partner(i) {}

  // Bond partner of ring-closing bond
  AtomIndex partner;
  // Marked number in depth-first traversal
  boost::optional<unsigned> number;
};

// Vertex data of atoms in the spanning tree
struct VertexProperties {
  unsigned rank = 0;
  bool aromatic = false;
  unsigned subsumedHydrogens = 0;
  std::vector<RingClosure> closures;
};

// Data type for representing the molecular graph broken down into a DAG
using SpanningTree = boost::adjacency_list<
  boost::vecS,
  boost::vecS,
  boost::bidirectionalS,
  VertexProperties
>;

namespace {

// Dummy weight map for use with Prim's minimum spanning tree algorithm
struct Unity {
  using key_type = typename SpanningTree::edge_descriptor;
  using value_type = double;
  using reference = double;
  using category = boost::readable_property_map_tag;
};

// Equally weight all edges in Prim's MST
inline double get(const Unity& /* u */, const Unity::key_type& /* e */) {
  return 1.0;
}

/* Weight map type for use with Prim's MST to bias customize spanning tree
 * generation towards SMILES purposes
 */
struct InverseBondOrder {
  using key_type = PrivateGraph::Edge;
  using value_type = double;
  using reference = double;
  using category = boost::readable_property_map_tag;

  explicit InverseBondOrder(const PrivateGraph& g) : graphRef(g) {}

  std::reference_wrapper<const PrivateGraph> graphRef;
};

/* Associated get function for InverseBondOrder, weighting:
 * - edges inversely proportional to the bond order
 * - heavy edge endpoint atoms proportional to their connectivity
 */
inline double get(const InverseBondOrder& u, const PrivateGraph::Edge& e) {
  const PrivateGraph& g = u.graphRef.get();
  const BondType type = g.bondType(e);
  const double order = Bond::bondOrderMap.at(static_cast<unsigned>(type));
  const unsigned inverseOrder = 7 - order;

  /* Weight edges by the graph degree (of heavy atoms only) to incentivize ring
   * closures at fulcrums of the heavy atom graph. Important for normalization
   * of e.g. biphenyl.
   */
  const auto heavyDegree = [&](const AtomIndex i) -> unsigned {
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

/* Returns references to static maps between permutation indices of
 * a vertex ordering for particular shapes to smiles stereodescriptor indices
 */
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

// Returns the smiles stereodescriptor matching a shape's vertex ordering
unsigned smilesStereodescriptor(
  const Shapes::Shape shape,
  const std::vector<Shapes::Vertex>& order
) {
  const auto allRotations = Shapes::Properties::generateAllRotations(shape, order);
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

// Helper struct performing smiles string construction
struct Emitter {
  //! Bond stereo markers in smiles
  enum class BondStereo { Forward, Backward };

  // SMILES has its own definition of what constitutes a heteroatom
  static bool heteroatom(const Utils::ElementType e) {
    return e != Utils::ElementType::H && e != Utils::ElementType::C;
  }

  Emitter(const Molecule& mol) : spanning(mol.V()), molecule(mol) {
    /* Plan:
     * - Find a minimum spanning tree while weighting edges according to their
     *   inverse bond order and sum of endpoint heavy atom degrees to find a
     *   suitable acyclization of the molecular graph for the smiles
     * - Find the longest path within that acyclic graph
     * - Find a minimum spanning tree starting at the root of the longest path
     *   and reduce edges from the (bidirectional acyclic) graph to fit that
     *   tree
     *
     * Details:
     * - Subsumed hydrogen atoms are disconnected in the construction of the
     *   initial bidirectional graph mimicking the molecule
     * - Record ring closures after acyclization
     */
    auto subsumption = markSubsumedHydrogens();
    const PrivateGraph& inner = mol.graph().inner();

    /* TODO maybe should place the root vertex of this MST at the 'centroid'
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

    // Add edges from predecessor list, but bidirectional!
    for(const PrivateGraph::Vertex v : inner.vertices()) {
      const PrivateGraph::Vertex predecessor = predecessors.at(v);

      // Skip unreachable vertices and the root
      if(predecessor == v) {
        continue;
      }

      // Skip subsumed hydrogen atoms
      if(subsumption.count(v) > 0) {
        continue;
      }

      boost::add_edge(predecessor, v, spanning);
      boost::add_edge(v, predecessor, spanning);
    }

    /* Figure out ring-closing bonds from predecessor list: If neither of the
     * bonds' vertices has the other atom as its predecessor, it is not in the
     * minimum spanning tree and represents a ring closure
     */
    for(const auto& bond : inner.edges()) {
      const PrivateGraph::Vertex u = inner.source(bond);
      const PrivateGraph::Vertex v = inner.target(bond);
      if(predecessors.at(u) != v && predecessors.at(v) != u) {
        spanning[u].closures.emplace_back(v);
        spanning[v].closures.emplace_back(u);
      }
    }

    longestPath = findMainBranch();
    const AtomIndex root = longestPath.front();

    /* Now for the second minimum spanning tree calculation to figure out the
     * directionality of the tree.
     */
    boost::prim_minimum_spanning_tree(
      spanning,
      boost::make_iterator_property_map(
        predecessors.begin(),
        boost::get(boost::vertex_index, inner.bgl())
      ),
      boost::root_vertex(root).
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

    recordVertexRank(root);
    markAromaticAtoms();
    markBondStereo();
  }

  /* Finds and marks hydrogen atoms that can be subsumed into heavy atom labels
   *
   * Returns a map from hydrogen atoms to their heavy atom to be subsumed into.
   * Increments the heavy atom vertex property counting the number of subsumed
   * hydrogen atoms.
   */
  std::unordered_map<AtomIndex, AtomIndex> markSubsumedHydrogens() {
    std::unordered_map<AtomIndex, AtomIndex> result;

    const PrivateGraph& inner = molecule.graph().inner();

    /* Find terminal hydrogen atoms that can be rewritten into heavy atom
     * labels, implicit or explicit
     */
    for(const AtomIndex i : inner.vertices()) {
      if(
        inner.elementType(i) != Utils::ElementType::H
        || inner.degree(i) != 1
      ) {
        continue;
      }

      for(const AtomIndex j : inner.adjacents(i)) {
        // Cannot subsume the terminal hydrogen in another hydrogen atom
        if(inner.elementType(j) == Utils::ElementType::H) {
          continue;
        }

        result.emplace(i, j);
        ++spanning[j].subsumedHydrogens;
      }
    }

    return result;
  }

  // Figure out the longest path in the acyclic graph
  std::vector<AtomIndex> findMainBranch() const {
    const unsigned V = molecule.V();
    AtomIndex suggestedRoot = guessRoot();
    std::vector<AtomIndex> predecessors(V, suggestedRoot);
    std::vector<unsigned> distances(V, 0);
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
    auto path = PredecessorMap {predecessors}.path(partner);

    // Normalization: Begin on a heteroatom, if possible
    if(
      heteroatom(molecule.graph().elementType(partner))
      && !heteroatom(molecule.graph().elementType(suggestedRoot))
    ) {
      std::swap(suggestedRoot, partner);
      std::reverse(std::begin(path), std::end(path));
    }

    return path;
  }

  /* Record the number of descendants of each vertex in the spanning tree
   * as the internal rank property (we'll use this in descent ordering later
   * for more normalized smiles).
   */
  void recordVertexRank(AtomIndex i) {
    // Do a non-tail-recursive DFS-like traversal
    for(
      auto iters = boost::out_edges(i, spanning);
      iters.first != iters.second;
      ++iters.first
    ) {
      const AtomIndex child = boost::target(*iters.first, spanning);
      recordVertexRank(child);
      spanning[i].rank += 1 + spanning[child].rank;
    }
  }

  // Set aromatic vertex properties of aromatic atoms
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
        spanning[i].aromatic = true;
      }
      for(const BondIndex& b : cycleEdges) {
        aromaticBonds.insert(b);
      }
    }
  }

  // Sets smiles descriptors on graph vertices for E/Z stereocenters
  void markBondStereo(const BondStereopermutator& permutator) {
    const BondIndex& bond = permutator.placement();

    if(
      permutator.numAssignments() < 2
      || molecule.bondType(bond) != BondType::Double
      || aromaticBonds.count(bond) > 0
    ) {
      return;
    }

    const auto orderedPlacement = [&]() {
      const auto forwardEdge = boost::edge(bond.first, bond.second, spanning);
      if(forwardEdge.second) {
        return std::make_pair(bond.first, bond.second);
      }
      return std::make_pair(bond.second, bond.first);
    }();

    const AtomIndex front = parent(orderedPlacement.first);

    if(molecule.bondType(BondIndex {front, orderedPlacement.first}) != BondType::Single) {
      return;
    }

    auto pathIter = std::find(
      std::begin(longestPath),
      std::end(longestPath),
      orderedPlacement.second
    );
    if(pathIter == std::end(longestPath)) {
      pathIter = std::begin(longestPath);
    }
    const AtomIndex back = descendants(orderedPlacement.second, pathIter).back();

    if(molecule.bondType(BondIndex {orderedPlacement.second, back}) != BondType::Single) {
      return;
    }

    const auto permutators = Temple::map(
      orderedPlacement,
      [&](const AtomIndex i) {
        return molecule.stereopermutators().option(i);
      }
    );

    const auto allowedShape = [](const Shapes::Shape shape) -> bool {
      return (
        shape == Shapes::Shape::Bent
        || shape == Shapes::Shape::EquilateralTriangle
      );
    };

    if(
      !permutators.first || !permutators.second
      || !allowedShape(permutators.first->getShape())
      || !allowedShape(permutators.second->getShape())
    ) {
      return;
    }

    const double dihedral = permutator.dihedral(
      permutators.first.value(),
      permutators.first->getRanking().getSiteIndexOf(front),
      permutators.second.value(),
      permutators.second->getRanking().getSiteIndexOf(back)
    );


    bondStereoMarkers.emplace(orderedPlacement.first, BondStereo::Forward);
    if(std::fabs(dihedral) < 1e-2) {
      // Front and back are Z to one another
      bondStereoMarkers.emplace(back, BondStereo::Backward);
    } else {
      // Front and back are E to one another
      bondStereoMarkers.emplace(back, BondStereo::Forward);
    }
  }

  void markBondStereo() {
    for(
      const BondStereopermutator& permutator :
      molecule.stereopermutators().bondStereopermutators()
    ) {
      markBondStereo(permutator);
    }
  }

  // Returns a flat map of distances to the argument vertex of the spanning tree
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

  // Guesses a possibly useful root of the MST by maximizing graph distance
  AtomIndex guessRoot() const {
    using IndexDistancePair = std::pair<AtomIndex, unsigned>;
    return Temple::accumulate(
      molecule.graph().atoms(),
      IndexDistancePair {0, 0},
      [&](const IndexDistancePair& carry, const AtomIndex i) -> IndexDistancePair {
        /* NOTE: We're using distances in the acyclic graph here, not distances
         * in the (possibly cyclic) molecule
         */
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

  // Construct the smiles representation of an atom
  std::string atomSymbol(
    const AtomIndex v,
    const std::vector<AtomIndex>& descendantOrder
  ) const {
    /* Decide whether to emit either
     *
     * - an aliphatic organic symbol
     * - an aromatic organic symbol
     * - or a bracket atom
     *   - This is necessary if isotopic or non-standard hydrogen count
     */
    const Utils::ElementType e = molecule.elementType(v);
    const bool isIsotope = Utils::ElementInfo::base(e) != e;
    const unsigned hydrogenCount = spanning[v].subsumedHydrogens;

    const bool isStereogenic = isSmilesStereogenic(v);

    if(!isIsotope && !isStereogenic && isValenceFillElement(e)) {
      const bool isAromatic = spanning[v].aromatic;

      const int valence = Temple::accumulate(
        molecule.graph().adjacents(v),
        0,
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

      // TODO need to kekulize here too, this doesn't quite work!

      const int aromaticityCompensatedValence = isAromatic ? valence + 1 : valence;

      const unsigned valenceFill = valenceFillElementImplicitHydrogenCount(
        aromaticityCompensatedValence,
        e
      );

      if(valenceFill == hydrogenCount) {
        // We can emit an aliphatic or aromatic organic symbol
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

    // Otherwise, bracket atom it is
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
    auto permutator = molecule.stereopermutators().option(v);
    assert(permutator);
    /* - Use the shape map to place the substituent order index of the smiles
     *   at the shape vertices (via the site indices), taking extra care of
     *   hydrogen counts (which are ordered first in the smiles sequence)
     * - Generate all rotations of that permutation and find a smiles
     *   stereodescriptor for that shape that matches
     */
    const auto& ranking = permutator->getRanking();
    std::vector<SiteIndex> siteSequence;
    if(spanning[v].subsumedHydrogens > 0) {
      // Subsumed hydrogens are placed first in the ordering (SMILES rules)
      SiteIndex i {0};
      for(const auto& siteAtoms : ranking.sites) {
        if(siteAtoms.size() == 1) {
          const AtomIndex siteAtom = siteAtoms.front();
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
    // Parent node
    for(auto iters = boost::in_edges(v, spanning); iters.first != iters.second; ++iters.first) {
      const AtomIndex parent = boost::source(*iters.first, spanning);
      siteSequence.push_back(ranking.getSiteIndexOf(parent));
    }
    // Ring closures (part of string after atom symbol)
    for(const RingClosure& closure : spanning[v].closures) {
      siteSequence.push_back(ranking.getSiteIndexOf(closure.partner));
    }
    // Descendants (explored later in DFS)
    for(const AtomIndex i : descendantOrder) {
      siteSequence.push_back(ranking.getSiteIndexOf(i));
    }

    assert(siteSequence.size() == Shapes::size(permutator->getShape()));

    const auto orderAtShapeVertices = Temple::map(
      siteSequence,
      Temple::Functor::at(permutator->getShapePositionMap())
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
      default: {
        std::cerr << "Precondition violation for smiles emitter atom stereodescriptor! Exiting.\n";
        std::exit(1); // Unreachable by precondition of isStereogenic
      }
    }
  }

  void writeEdgeType(const AtomIndex v, const AtomIndex w) {
    const BondIndex bond {v, w};
    switch(molecule.bondType(bond)) {
      case BondType::Single: {
        /* Bond between aromatic atoms that isn't aromatic itself,
         * e.g. in biphenyl c1ccccc1-c2ccccc2
         */
        if(
          spanning[v].aromatic
          && spanning[w].aromatic
          && aromaticBonds.count(bond) == 0
        ) {
          smiles += "-";
        }

        /* Double bond stereo markers */
        const auto markerIter = bondStereoMarkers.find(w);
        if(markerIter != std::end(bondStereoMarkers)) {
          switch(markerIter->second) {
            case BondStereo::Forward: {
              smiles += "/";
              break;
            }
            case BondStereo::Backward: {
              smiles += "\\";
              break;
            }
          }
        }

        break;
      }
      case BondType::Double: { smiles += "="; break; }
      case BondType::Triple: { smiles += "#"; break; }
      // NOTE: Eta and bond orders four through six aren't explicitly written
      default: break;
    }
  }

  /* Sort descendant vertices by increasing number of descendants, ensuring
   * the continuation of the longest path is always last.
   *
   * TODO check if the additional conditions regarding the main branch are
   * actually necessary. In principle, it could be guaranteed that the longest
   * path is last by virtue of being the longest path (i.e. the path with the
   * most descendants). It might not be stable against multiple longest paths,
   * though.
   */
  std::vector<AtomIndex> descendants(
    const AtomIndex v,
    const std::vector<AtomIndex>::iterator& pathIter
  ) const {
    const auto nextOnMainBranch = [&]() -> boost::optional<AtomIndex> {
      if(v != *pathIter) {
        return boost::none;
      }

      const auto nextIter = pathIter + 1;
      if(nextIter == std::end(longestPath)) {
        return boost::none;
      }

      return *nextIter;
    }();

    std::vector<AtomIndex> vertices;
    for(auto iters = boost::out_edges(v, spanning); iters.first != iters.second; ++iters.first) {
      vertices.push_back(boost::target(*iters.first, spanning));
    }

    Temple::sort(
      vertices,
      [&](const AtomIndex i, const AtomIndex j) -> bool {
        if(i == nextOnMainBranch) {
          return false;
        }

        if(j == nextOnMainBranch) {
          return true;
        }

        return spanning[i].rank < spanning[j].rank;
      }
    );

    return vertices;
  }

  void dfs(const AtomIndex v, std::vector<AtomIndex>::iterator& pathIter) {
    const auto descendantOrder = descendants(v, pathIter);
    smiles += atomSymbol(v, descendantOrder);
    // Assign ring closing bonds numbers
    for(RingClosure& closure : spanning[v].closures) {
      if(!closure.number) {
        closure.number = closureIndex;

        // Find matching closure for partner
        auto& partnerClosures = spanning[closure.partner].closures;
        const auto partnerIter = std::find_if(
          std::begin(partnerClosures),
          std::end(partnerClosures),
          [&](const RingClosure& hay) {
            return hay.partner == v;
          }
        );
        assert(partnerIter != std::end(partnerClosures));
        assert(!partnerIter->number);

        // Mark partner's number too
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

    // Descend into substituents
    for(const AtomIndex w : descendantOrder) {
      const bool last = (w == descendantOrder.back());
      if(!last) {
        smiles += "(";
      }
      writeEdgeType(v, w);
      dfs(w, pathIter);
      if(!last) {
        smiles += ")";
      }
    }
  }

  AtomIndex parent(const AtomIndex v) const {
    auto iters = boost::in_edges(v, spanning);
    assert(iters.first != iters.second);
    return boost::source(*iters.first, spanning);
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
  /* NOTE: This cannot be a spanning-tree edge property since aromatic ring
   * closures are not edges in the tree
   */
  std::unordered_set<BondIndex, boost::hash<BondIndex>> aromaticBonds;
  std::unordered_map<AtomIndex, BondStereo> bondStereoMarkers;
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
