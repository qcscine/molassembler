/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Graph/EditDistance.h"

// #include "boost/graph/graphviz.hpp"

#include <Eigen/Dense>
#include <queue>
#include <deque>

namespace Scine {
namespace Molassembler {
namespace GraphAlgorithms {

constexpr PrivateGraph::Vertex EditDistanceForest::epsilon;

double EditDistanceForest::HausdorffHeuristic::edgeFunction(
  const PrivateGraph::Vertex u,
  const PrivateGraph::Vertex v
) const {
  double cost = 0;
  for(const PrivateGraph::Edge aEdge : a.edges(u)) {
    double minEdgeCost = costFn.edgeAlteration();
    for(PrivateGraph::Edge bEdge : b.edges(v)) {
      minEdgeCost = std::min(
        minEdgeCost,
        costFn.bondSubstitution(a.bondType(aEdge), b.bondType(bEdge)) / 2.0
      );
    }
    cost += minEdgeCost;
  }

  for(const PrivateGraph::Edge bEdge : b.edges(v)) {
    double minEdgeCost = costFn.edgeAlteration();
    for(PrivateGraph::Edge aEdge : a.edges(u)) {
      minEdgeCost = std::min(
        minEdgeCost,
        costFn.bondSubstitution(a.bondType(aEdge), b.bondType(bEdge)) / 2.0
      );
    }
    cost += minEdgeCost;
  }

  return std::max(
    cost,
    std::fabs(a.degree(u) - b.degree(v)) * costFn.edgeAlteration()
  );
}

double EditDistanceForest::HausdorffHeuristic::nodeFunction(
  const PrivateGraph::Vertex u,
  const PrivateGraph::Vertex v
) const {
  if(u == epsilon) {
    return costFn.vertexAlteration() + b.degree(v) * costFn.edgeAlteration() / 2.0;
  }

  if(v == epsilon) {
    return costFn.vertexAlteration() + a.degree(u) * costFn.edgeAlteration() / 2.0;
  }

  return (
    costFn.elementSubstitution(a.elementType(u), b.elementType(v))
    + edgeFunction(u, v) / 2.0
  ) / 2.0;
}

double EditDistanceForest::HausdorffHeuristic::estimate(
  const unsigned depth,
  const std::unordered_set<PrivateGraph::Vertex>& usedBVertices
) const {
  // Vertices free in a are those greater or equal to depth
  // Vertices free in b are those not in usedBVertices
  std::unordered_map<PrivateGraph::Vertex, double> aCost;
  const PrivateGraph::Vertex aN = a.V();
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

EditDistanceForest::Vertex
EditDistanceForest::search(const PrivateGraph& a, const PrivateGraph& b) {
  const unsigned completeDepth = std::max(a.V(), b.V());

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

    if(traversal.depth < a.V()) {
      std::unordered_set<unsigned> matchedComponents;
      const auto addBranch = [&](const PrivateGraph::Vertex v) {
        if(bVertices.count(v) > 0) {
          return;
        }
        const unsigned component = bHydrogenComponents.at(v);
        if(matchedComponents.count(component) > 0) {
          return;
        }
        matchedComponents.insert(component);

        const Utils::ElementType ukElementType = a.elementType(traversal.depth);
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
          if(aBondOption && bBondOption) {
            impliedEdgeCosts += costFn.bondSubstitution(*aBondOption, *bBondOption);
          } else if(static_cast<bool>(aBondOption) ^ static_cast<bool>(bBondOption)) {
            impliedEdgeCosts += costFn.edgeAlteration();
          }

          ++i;
        }
        const unsigned cost = (
          g[top].costSum
          + costFn.elementSubstitution(ukElementType, b.elementType(v))
          + impliedEdgeCosts
        );
        const Vertex newVertex = addVertex(
          VertexData {
            v,
            cost,
            cost + estimate(traversal, bVertices, v)
          },
          top
        );
        queue.push(newVertex);
      };

      const auto preconditionIter = preconditioning.left.find(traversal.depth);
      if(preconditionIter != preconditioning.left.end()) {
        addBranch(preconditionIter->second);
      } else {
        for(const PrivateGraph::Vertex v : b.vertices()) {
          addBranch(v);
        }
        // Vertex deletion
        unsigned impliedEdgeCosts = 0;
        for(PrivateGraph::Vertex i = 0; i < traversal.depth; ++i) {
          const auto aEdgeOption = a.edgeOption(i, traversal.depth);
          if(aEdgeOption) {
            // If the implied edge exists, it has to be deleted
            impliedEdgeCosts += costFn.edgeAlteration();
          }
        }
        const unsigned cost = g[top].costSum + costFn.vertexAlteration() + impliedEdgeCosts;
        const Vertex newVertex = addVertex(
          VertexData {
            epsilon,
            cost,
            cost + estimate(traversal, bVertices, epsilon)
          },
          top
        );
        queue.push(newVertex);
      }
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
            impliedEdgeCosts += costFn.edgeAlteration();
          }
        }
        const unsigned cost = g[currentForestVertex].costSum + costFn.vertexAlteration() + impliedEdgeCosts;
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
  }
}

} // namespace GraphAlgorithms
} // namespace Molassembler
} // namespace Scine
