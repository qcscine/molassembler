/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Graph/EditDistance.h"

// #include "boost/graph/graphviz.hpp"

#include <Eigen/Dense>
#include <limits>
#include <queue>
#include <deque>

namespace Scine {
namespace Molassembler {
namespace GraphAlgorithms {
namespace {

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

  /* Cost function summary
   * Vertex labeling cost = 1 if element type mismatch
   * Edge labeling cost = 1 if bond type mismatch (maaaybe distance metric too?)
   *
   * Vertex insert / delete = 1
   */
  struct Cost {
    static constexpr unsigned insertion = 1;
    static constexpr unsigned deletion = 1;

    static constexpr unsigned substitution(Utils::ElementType from, Utils::ElementType to) {
      if(from != to) {
        return 1;
      }

      return 0;
    }

    static constexpr unsigned substitution(BondType from, BondType to) {
      if(from != to) {
        return 1;
      }

      return 0;
    }
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
    double edgeFunction(
      const PrivateGraph::Vertex u,
      const PrivateGraph::Vertex v
    ) const {
      double cost = 0;
      for(const PrivateGraph::Edge aEdge : a.edges(u)) {
        double minEdgeCost = Cost::deletion;
        for(PrivateGraph::Edge bEdge : b.edges(v)) {
          minEdgeCost = std::min(
            minEdgeCost,
            Cost::substitution(a.bondType(aEdge), b.bondType(bEdge)) / 2.0
          );
        }
        cost += minEdgeCost;
      }

      for(const PrivateGraph::Edge bEdge : b.edges(v)) {
        double minEdgeCost = Cost::deletion;
        for(PrivateGraph::Edge aEdge : a.edges(u)) {
          minEdgeCost = std::min(
            minEdgeCost,
            Cost::substitution(a.bondType(aEdge), b.bondType(bEdge)) / 2.0
          );
        }
        cost += minEdgeCost;
      }

      return std::max(
        cost,
        std::fabs(a.degree(u) - b.degree(v)) * Cost::deletion
      );
    }

    double nodeFunction(
      const PrivateGraph::Vertex u,
      const PrivateGraph::Vertex v
    ) const {
      if(u == epsilon) {
        return Cost::insertion + b.degree(v) * Cost::insertion / 2.0;
      }

      if(v == epsilon) {
        return Cost::deletion + a.degree(u) * Cost::deletion / 2.0;
      }

      return (
        Cost::substitution(a.elementType(u), b.elementType(v))
        + edgeFunction(u, v) / 2.0
      ) / 2.0;
    }

    HausdorffHeuristic(const PrivateGraph& i, const PrivateGraph& j) : a(i), b(j) {
      precomputedNodeFunction.resize(a.N(), b.N());
      for(PrivateGraph::Vertex u : a.vertices()) {
        for(PrivateGraph::Vertex v : b.vertices()) {
          precomputedNodeFunction(u, v) = nodeFunction(u, v);
        }
      }
    }

    double estimate(
      const unsigned depth,
      const std::unordered_set<PrivateGraph::Vertex>& usedBVertices
    ) const {
      // Vertices free in a are those greater or equal to depth
      // Vertices free in b are those not in usedBVertices
      std::unordered_map<PrivateGraph::Vertex, double> aCost;
      const PrivateGraph::Vertex aN = a.N();
      for(PrivateGraph::Vertex u = depth; u < aN; ++u) {
        aCost.emplace(u, nodeFunction(u, epsilon));
      }

      std::unordered_map<PrivateGraph::Vertex, double> bCost;
      for(PrivateGraph::Vertex v : b.vertices()) {
        if(usedBVertices.count(v) > 0) {
          continue;
        }
        bCost.emplace(v, nodeFunction(epsilon, v));
      }

      for(auto& aCostPair : aCost) {
        for(auto& bCostPair : bCost) {
          const double cost = precomputedNodeFunction(aCostPair.first, bCostPair.first);
          aCostPair.second = std::min(aCostPair.second, cost);
          bCostPair.second = std::min(bCostPair.second, cost);
        }
      }

      double sum = 0.0;
      for(auto aCostPair : aCost) {
        sum += aCostPair.second;
      }
      for(auto bCostPair : bCost) {
        sum += bCostPair.second;
      }
      return sum;
    }

    const PrivateGraph& a;
    const PrivateGraph& b;
    Eigen::MatrixXd precomputedNodeFunction;
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
//!@}

  EditDistanceForest(const PrivateGraph& a, const PrivateGraph& b)
    : queue(QueueComparator(g)), heuristic(a, b)
  {
    // Initialize the search forest
    std::vector<Vertex> roots;
    if(a.N() > 0) {
      const Utils::ElementType u1ElementType = a.elementType(0);
      for(const PrivateGraph::Vertex v : b.vertices()) {
        Vertex newVertex = boost::add_vertex(g);
        const unsigned cost = Cost::substitution(u1ElementType, b.elementType(v));
        g[newVertex] = VertexData {
          v,
          cost,
          cost + heuristic.estimate(1, {v})
        };
        roots.push_back(newVertex);
      }
      Vertex deletion = boost::add_vertex(g);
      const unsigned cost = Cost::deletion;
      g[deletion] = VertexData {
        epsilon,
        cost,
        cost + heuristic.estimate(1, {})
      };
      roots.push_back(deletion);
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
   * @brief A* search algorithm for minimal GED
   *
   * Kaspar Riesen, Structural Pattern Recognition with Graph Edit
   * Distance, (Springer) 2015, pp. 35ff,
   * https://doi.org/10.1007/978-3-319-27252-8
   */
  Vertex search(const PrivateGraph& a, const PrivateGraph& b) {
    const unsigned completeDepth = std::max(a.N(), b.N());

    // unsigned counter = 0;

    // A* algorithm
    while(true) {
      // Select the current best path
      const Vertex top = queue.top();
      queue.pop();

      const Traversal traversal = traverse(top);
      std::unordered_set<PrivateGraph::Vertex> bVertices {
        std::begin(traversal.bVertices),
        std::end(traversal.bVertices)
      };
      bVertices.erase(epsilon);
      if(traversal.depth == completeDepth) {
        // Current node is a complete edit path
        return top;
      }

      if(traversal.depth < a.N()) {
        for(const PrivateGraph::Vertex v : b.vertices()) {
          if(bVertices.count(v) > 0) {
            continue;
          }

          const Utils::ElementType ukElementType = a.elementType(traversal.depth);
          const unsigned vertexCost = ukElementType == b.elementType(v) ? 0U : 1U;
          /* Calculate edge costs of all implied edges
           * - New vertex mapping between graphs in this loop is depth -> v
           * - If we loop over all pairings in this branch i -> *rIter,
           *   then all implied edges are (i-depth) -> (*rIter-v)
           */
          unsigned impliedEdgeCosts = 0;
          PrivateGraph::Vertex i = 0;
          const auto rEnd = std::rend(traversal.bVertices);
          for(auto rIter = std::rbegin(traversal.bVertices); rIter != rEnd; ++rIter) {
            const auto aBondOption = bondTypeOption(i, traversal.depth, a);
            const auto bBondOption = bondTypeOption(*rIter, v, b);
            if(aBondOption != bBondOption) {
              impliedEdgeCosts += 1;
            }

            ++i;
          }
          const unsigned cost = g[top].costSum + vertexCost + impliedEdgeCosts;
          const Vertex newVertex = addVertex(
            VertexData {
              v,
              cost,
              cost + estimate(traversal, bVertices, v)
            },
            top
          );
          queue.push(newVertex);
        }
        // Vertex deletion
        unsigned impliedEdgeCosts = 0;
        for(PrivateGraph::Vertex i = 0; i < traversal.depth; ++i) {
          const auto aEdgeOption = a.edgeOption(i, traversal.depth);
          if(aEdgeOption) {
            // If the implied edge exists, it has to be deleted
            impliedEdgeCosts += Cost::deletion;
          }
        }
        const unsigned cost = g[top].costSum + Cost::deletion + impliedEdgeCosts;
        const Vertex newVertex = addVertex(
          VertexData {
            epsilon,
            cost,
            cost + estimate(traversal, bVertices, epsilon)
          },
          top
        );
        queue.push(newVertex);
      } else {
        // Add all vertex insertions at once
        Vertex currentForestVertex = top;
        for(PrivateGraph::Vertex v : b.vertices()) {
          if(bVertices.count(v) > 0) {
            continue;
          }

          unsigned impliedEdgeCosts = 0;
          for(const PrivateGraph::Vertex w : bVertices) {
            const auto bEdgeOption = b.edgeOption(v, w);
            if(bEdgeOption) {
              impliedEdgeCosts += Cost::insertion;
            }
          }
          const unsigned cost = g[currentForestVertex].costSum + Cost::insertion + impliedEdgeCosts;
          /* NOTE: No point in estimating here since the final value will be
           * true and no intermediate solutions will be added to the queue, so
           * we skip it
           */
          currentForestVertex = addVertex(
            VertexData {
              v,
              cost,
              static_cast<double>(cost)
            },
            currentForestVertex
          );
          bVertices.insert(v);
        }
        queue.push(currentForestVertex);
      }

      // std::string fname = "step-" + std::to_string(counter) + ".dot";
      // std::ofstream fout(fname);
      // boost::write_graphviz(fout, g, GraphvizWriter {g});
      // ++counter;
    }
  }

//!@name Data members
//!@{
  Forest g;
  Queue queue;
  HausdorffHeuristic heuristic;
  Vertex result;
//!@}
};

constexpr PrivateGraph::Vertex EditDistanceForest::epsilon;

} // namespace

unsigned editDistance(const PrivateGraph& a, const PrivateGraph& b) {
  EditDistanceForest forest {a, b};
  return forest.g[forest.result].costSum;
}

} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine
