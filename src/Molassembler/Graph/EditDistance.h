/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#ifndef INCLUDE_MOLASSEMBLER_GRAPH_EDIT_DISTANCE_H
#define INCLUDE_MOLASSEMBLER_GRAPH_EDIT_DISTANCE_H

#include "Molassembler/GraphAlgorithms.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/Subgraphs.h"

#include <Eigen/Dense>
#include <queue>
#include <deque>

namespace Scine {
namespace Molassembler {
namespace GraphAlgorithms {

/**
 * @brief Forest data structure to calculate exact graph edit distance with
 */
struct EditDistanceForest {
//!@name Types
//!@{
  //! Data stored for each vertex in a tree
  struct VertexData {
    PrivateGraph::Vertex mappedVertex;
    unsigned costSum;
    double estimatedCost;
  };

  //! Deleted vertex sentinel value
  static constexpr PrivateGraph::Vertex epsilon = std::numeric_limits<AtomIndex>::max();

  using Forest = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::bidirectionalS,
    VertexData
  >;
  using Vertex = Forest::vertex_descriptor;

  struct QueueComparator {
    explicit QueueComparator(const Forest& x) : forestPtr(&x) {}

    bool operator() (const Vertex a, const Vertex b) const {
      const Forest& t = *forestPtr;
      /* Note: priority queue outputs 'largest' items first, so we want items
       * with the lowest cost sum last. This means using the inverse comparison
       * as below:
       */
      return t[a].estimatedCost > t[b].estimatedCost;
    }

    const Forest* forestPtr;
  };

  using Queue = std::priority_queue<Vertex, std::deque<Vertex>, QueueComparator>;

  struct GraphvizWriter {
    explicit GraphvizWriter(const Forest& x) : forest(x) {}

    void operator() (std::ostream& os, Vertex v) const {
      const VertexData& data = forest[v];
      os << "[label=\"";
      if(data.mappedVertex == epsilon) {
        os << "eps";
      } else {
        os << data.mappedVertex;
      }
      os << ", (" << data.estimatedCost <<  ", " << data.costSum << ")\"]";
    }

    const Forest& forest;
  };

  struct Traversal {
    unsigned depth;
    std::vector<PrivateGraph::Vertex> bVertices;
  };

  /**
   * @brief Heuristic to estimate future cost of node matching
   *
   * Fischer A., Plamondon R., Savaria Y., Riesen K., Bunke H. (2014) A
   * Hausdorff Heuristic for Efficient Computation of Graph Edit Distance. In:
   * Fränti P., Brown G., Loog M., Escolano F., Pelillo M. (eds) Structural,
   * Syntactic, and Statistical Pattern Recognition. S+SSPR 2014. Lecture Notes
   * in Computer Science, vol 8621. Springer, Berlin, Heidelberg.
   * https://doi.org/10.1007/978-3-662-44415-3_9
   */
  struct HausdorffHeuristic {
    HausdorffHeuristic(
      const PrivateGraph& i,
      const PrivateGraph& j,
      const EditCost& cost
    ) : a(i), b(j), costFn(cost) {
      precomputedNodeFunction.resize(a.V(), b.V());
      for(PrivateGraph::Vertex u : a.vertices()) {
        for(PrivateGraph::Vertex v : b.vertices()) {
          precomputedNodeFunction(u, v) = nodeFunction(u, v);
        }
      }
    }

    double edgeFunction(PrivateGraph::Vertex u, PrivateGraph::Vertex v) const;
    double nodeFunction(PrivateGraph::Vertex u, PrivateGraph::Vertex v) const;

    double estimate(
      unsigned depth,
      const std::unordered_set<PrivateGraph::Vertex>& usedBVertices
    ) const;

    const PrivateGraph& a;
    const PrivateGraph& b;
    Eigen::MatrixXd precomputedNodeFunction;
    const EditCost& costFn;
  };
//!@}

//!@name Static methods
//!@{
  static boost::optional<BondType> bondTypeOption(
    const PrivateGraph::Vertex a,
    const PrivateGraph::Vertex b,
    const PrivateGraph& g
  ) {
    if(a == epsilon || b == epsilon) {
      return boost::none;
    }

    if(auto edgeOption = g.edgeOption(a, b)) {
      return g.bondType(edgeOption.value());
    }

    return boost::none;
  }

  static std::vector<unsigned> hydrogenComponents(const PrivateGraph& graph) {
    auto isTerminalHydrogen = [&](AtomIndex v) {
      return graph.elementType(v) == Utils::ElementType::H && graph.degree(v) == 1;
    };

    std::vector<unsigned> components(graph.V());
    unsigned component = 0;
    for(const AtomIndex i : graph.vertices()) {
      if(isTerminalHydrogen(i)) {
        continue;
      }

      components.at(i) = component;
      ++component;

      bool anyHydrogensMarked = false;
      for(const AtomIndex j : graph.adjacents(i)) {
        if(isTerminalHydrogen(j)) {
          components.at(j) = component;
          anyHydrogensMarked = true;
        }
      }
      if(anyHydrogensMarked) {
        ++component;
      }
    }

    return components;
  }
//!@}
  EditDistanceForest(
    const PrivateGraph& a,
    const PrivateGraph& b,
    const EditCost& costFnRef,
    const Subgraphs::IndexMap& preconditioningMap
  ) : preconditioning(preconditioningMap),
      bHydrogenComponents(hydrogenComponents(b)),
      queue(QueueComparator(g)),
      costFn(costFnRef),
      heuristic(a, b, costFn)
  {
    // Initialize the search forest
    std::vector<Vertex> roots;

    if(a.V() > 0) {
      const Utils::ElementType u1ElementType = a.elementType(0);
      std::unordered_set<unsigned> matchedComponents;
      const auto addRoot = [&](const PrivateGraph::Vertex v) {
        const unsigned component = bHydrogenComponents.at(v);
        if(matchedComponents.count(component) > 0) {
          return;
        }
        matchedComponents.insert(component);

        Vertex newVertex = boost::add_vertex(g);
        const unsigned cost = costFn.elementSubstitution(u1ElementType, b.elementType(v));
        g[newVertex] = VertexData {
          v,
          cost,
          cost + heuristic.estimate(1, {v})
        };
        roots.push_back(newVertex);
      };

      const auto preconditionIter = preconditioning.left.find(0);
      if(preconditionIter != preconditioning.left.end()) {
        // First vertex is preconditioned. Only add it to roots and nothing else
        addRoot(preconditionIter->second);
      } else {
        // No preconditions, handle all possible matches and deletion
        for(const PrivateGraph::Vertex v : b.vertices()) {
          addRoot(v);
        }
        Vertex deletion = boost::add_vertex(g);
        const unsigned cost = costFn.vertexAlteration();
        g[deletion] = VertexData {
          epsilon,
          cost,
          cost + heuristic.estimate(1, {})
        };
        roots.push_back(deletion);
      }
    }

    queue = Queue {std::begin(roots), std::end(roots), QueueComparator(g)};
    result = search(a, b);
  }

  Traversal traverse(Vertex v) const {
    Traversal traversal {1, {}};
    traversal.bVertices.push_back(g[v].mappedVertex);
    while(boost::in_degree(v, g) > 0) {
      auto inEdgesIterPair = boost::in_edges(v, g);
      assert(inEdgesIterPair.first != inEdgesIterPair.second);
      v = boost::source(*(inEdgesIterPair.first), g);
      traversal.depth += 1;
      traversal.bVertices.push_back(g[v].mappedVertex);
    }
    return traversal;
  };

  Vertex addVertex(VertexData data, Vertex parent) {
    Vertex newVertex = boost::add_vertex(g);
    g[newVertex] = std::move(data);
    boost::add_edge(parent, newVertex, g);
    return newVertex;
  }

  double estimate(
    const Traversal& traversal,
    std::unordered_set<PrivateGraph::Vertex>& bVertices,
    PrivateGraph::Vertex newTarget
  ) const {
    bVertices.insert(newTarget);
    const double value = heuristic.estimate(traversal.depth + 1, bVertices);
    bVertices.erase(newTarget);
    return value;
  }

  /**
   * @brief Best-first A* search algorithm for minimal GED
   *
   * Kaspar Riesen, Structural Pattern Recognition with Graph Edit
   * Distance, (Springer) 2015, pp. 35ff,
   * https://doi.org/10.1007/978-3-319-27252-8
   */
  Vertex search(const PrivateGraph& a, const PrivateGraph& b);

//!@name Data members
//!@{
  Forest g;
  const Subgraphs::IndexMap& preconditioning;
  std::vector<unsigned> bHydrogenComponents;
  Queue queue;
  const EditCost& costFn;
  HausdorffHeuristic heuristic;
  Vertex result;
//!@}
};

} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine

#endif
