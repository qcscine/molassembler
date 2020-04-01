/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Interpret.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"

#include "molassembler/BondOrders.h"
#include "molassembler/Molecule.h"
#include "molassembler/Graph.h"
#include "molassembler/Graph/PrivateGraph.h"

#include "molassembler/Temple/Functional.h"

namespace Scine {
namespace molassembler {
namespace interpret {

namespace {

Utils::PositionCollection paste(const std::vector<Utils::Position>& positions) {
  const unsigned N = positions.size();
  auto matrix = Utils::PositionCollection(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    matrix.row(i) = positions.at(i);
  }
  return matrix;
}

} // namespace

struct MoleculeParts {
  PrivateGraph graph;
  std::vector<Utils::Position> angstromPositions;
  boost::optional<
    std::vector<BondIndex>
  > bondStereopermutatorCandidatesOptional;
};

struct Parts {
  std::vector<MoleculeParts> precursors;
  std::vector<unsigned> componentMap;
  unsigned nZeroLengthPositions = 0;
};

// Yields a graph structure without element type annotations
PrivateGraph discretize(
  const Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization
) {
  const PrivateGraph::Vertex N = bondOrders.getSystemSize();
  PrivateGraph graph {N};

  if(discretization == BondDiscretizationOption::Binary) {
    for(unsigned i = 0; i < N; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        double bondOrder = bondOrders.getOrder(i, j);

        if(bondOrder > 0.5) {
          graph.addEdge(i, j, BondType::Single);
        }
      }
    }
  } else if(discretization == BondDiscretizationOption::RoundToNearest) {
    for(unsigned i = 0; i < N; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        double bondOrder = bondOrders.getOrder(i, j);

        if(bondOrder > 0.5) {
          auto bond = static_cast<BondType>(
            std::round(bondOrder) - 1
          );

          if(bondOrder > 6.5) {
            bond = BondType::Sextuple;
          }

          graph.addEdge(i, j, bond);
        }
      }
    }
  }

  return graph;
}

Parts construeParts(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromWrapper,
  const Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorThreshold
) {
  const unsigned N = elements.size();

  // Check preconditions
  if(angstromWrapper.positions.rows() != N) {
    throw std::invalid_argument(
      "Number of positions in angstrom wrapper do not match number of elements"
    );
  }

  if(bondOrders.getSystemSize<unsigned>() != N) {
    throw std::invalid_argument(
      "Bond order argument system size does not match number of elements"
    );
  }

  PrivateGraph atomCollectionGraph = discretize(bondOrders, discretization);

  Parts parts;
  const unsigned numComponents = atomCollectionGraph.connectedComponents(parts.componentMap);
  parts.precursors.resize(numComponents);

  if(stereopermutatorThreshold) {
    for(auto& precursor : parts.precursors) {
      // Empty-vector-initialize the candidate optionals
      precursor.bondStereopermutatorCandidatesOptional = std::vector<BondIndex> {};
    }
  }

  // Map from original index to component index
  std::vector<PrivateGraph::Vertex> indexInComponentMap (N);

  /* Maybe
   * - filtered_graph using predicate of componentMap number
   * - copy_graph to new Graph keeping element types and bond orders
   *
   * - alternately, must keep a map of atomcollection index to precursor index
   *   and new precursor atom index in order to transfer edges too
   */

  for(unsigned i = 0; i < N; ++i) {
    auto& precursor = parts.precursors.at(
      parts.componentMap.at(i)
    );

    // Add a new vertex with element information
    PrivateGraph::Vertex newIndex = precursor.graph.addVertex(elements.at(i));

    // Save new index in precursor graph
    indexInComponentMap.at(i) = newIndex;

    if(angstromWrapper.positions.row(i).norm() <= 1e-14) {
      parts.nZeroLengthPositions += 1;
    }

    // Copy over position information
    precursor.angstromPositions.emplace_back(
      angstromWrapper.positions.row(i)
    );
  }

  // Copy over edges and bond orders
  for(const PrivateGraph::Edge& edge : atomCollectionGraph.edges()) {
    const PrivateGraph::Vertex source = atomCollectionGraph.source(edge);
    const PrivateGraph::Vertex target = atomCollectionGraph.target(edge);

    // Both source and target are part of the same component (since they are bonded)
    auto& precursor = parts.precursors.at(
      parts.componentMap.at(source)
    );

    // Copy over the edge
    precursor.graph.addEdge(
      indexInComponentMap.at(source),
      indexInComponentMap.at(target),
      atomCollectionGraph.bondType(edge)
    );

    // If the edge's bond order exceeds the threshold optional
    if(
      stereopermutatorThreshold
      && bondOrders.getOrder(source, target) >= *stereopermutatorThreshold
    ) {
      precursor.bondStereopermutatorCandidatesOptional->emplace_back(
        indexInComponentMap.at(source),
        indexInComponentMap.at(target)
      );
    }
  }

  return parts;
}

MoleculesResult molecules(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromWrapper,
  const Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorThreshold
) {
  Parts parts = construeParts(
    elements,
    angstromWrapper,
    bondOrders,
    discretization,
    stereopermutatorThreshold
  );

  // Collect results
  MoleculesResult result;

  /* Transform precursors into Molecules. Positions may only be used if there
   * is at most one position very close to (0, 0, 0). Otherwise, we assume that
   * the given positions are faulty or no positional information is present,
   * and only the graph is used to create the Molecules.
   */
  if(parts.nZeroLengthPositions < 2) {
    result.molecules.reserve(parts.precursors.size());
    for(auto& precursor : parts.precursors) {
      result.molecules.emplace_back(
        Graph {std::move(precursor.graph)},
        AngstromPositions(paste(precursor.angstromPositions), LengthUnit::Angstrom),
        precursor.bondStereopermutatorCandidatesOptional
      );
    }
  } else {
    result.molecules.reserve(parts.precursors.size());
    for(auto& precursor : parts.precursors) {
      result.molecules.emplace_back(
        Graph {std::move(precursor.graph)}
      );
    }
  }

  result.componentMap = std::move(parts.componentMap);

  return result;
}

MoleculesResult molecules(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromWrapper,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorThreshold
) {
  return molecules(
    elements,
    angstromWrapper,
    uffBondOrders(elements, angstromWrapper),
    discretization,
    stereopermutatorThreshold
  );
}

MoleculesResult molecules(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorThreshold
) {
  return molecules(
    atomCollection.getElements(),
    AngstromPositions {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    discretization,
    stereopermutatorThreshold
  );
}

MoleculesResult molecules(
  const Utils::AtomCollection& atomCollection,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorThreshold
) {
  AngstromPositions angstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr};

  return molecules(
    atomCollection.getElements(),
    angstromWrapper,
    uffBondOrders(atomCollection.getElements(), angstromWrapper),
    discretization,
    stereopermutatorThreshold
  );
}

std::vector<Utils::AtomCollection> applyInterpretationMap(
  const ComponentMap& componentMap,
  const Utils::AtomCollection& atomCollection
) {
  const unsigned nComponents = *std::max_element(
    std::begin(componentMap),
    std::end(componentMap)
  ) + 1;
  std::vector<unsigned> componentSizes(nComponents, 0);
  for(unsigned i : componentMap) {
    componentSizes.at(i) += 1;
  }

  /* Allocate the collections */
  std::vector<Utils::AtomCollection> collections = temple::map(
    componentSizes,
    [](const unsigned size) { return Utils::AtomCollection(size); }
  );
  std::vector<unsigned> collectionSizeCount(nComponents, 0);

  for(unsigned i = 0; i < componentMap.size(); ++i) {
    unsigned moleculeIndex = componentMap.at(i);
    Utils::AtomCollection& collection = collections.at(moleculeIndex);
    unsigned& collectionSize = collectionSizeCount.at(moleculeIndex);
    collection.setElement(
      collectionSize,
      atomCollection.getElement(i)
    );
    collection.setPosition(
      collectionSize,
      atomCollection.getPosition(i)
    );
    ++collectionSize;
  }

  return collections;
}

std::vector<
  std::vector<unsigned>
> invertComponentMap(const ComponentMap& componentMap) {
  const unsigned nComponents = *std::max_element(
    std::begin(componentMap),
    std::end(componentMap)
  ) + 1;

  std::vector<
    std::vector<unsigned>
  > inverseMaps (nComponents);

  const unsigned N = componentMap.size();
  for(unsigned i = 0; i < N; ++i) {
    inverseMaps.at(componentMap.at(i)).push_back(i);
  }

  return inverseMaps;
}

GraphsResult graphs(
  const Utils::ElementTypeCollection& elements,
  const AngstromPositions& angstromWrapper,
  const Utils::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization
) {
  Parts parts = construeParts(
    elements,
    angstromWrapper,
    bondOrders,
    discretization,
    boost::none
  );

  GraphsResult result;
  result.graphs.reserve(parts.precursors.size());
  for(auto& precursor : parts.precursors) {
    result.graphs.emplace_back(std::move(precursor.graph));
  }
  result.componentMap = std::move(parts.componentMap);
  return result;
}

GraphsResult graphs(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders,
  BondDiscretizationOption discretization
) {
  return graphs(
    atomCollection.getElements(),
    AngstromPositions {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    discretization
  );
}

} // namespace interpret
} // namespace molassembler
} // namespace Scine
