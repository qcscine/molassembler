/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "molassembler/Interpret.h"

#include "boost/range/iterator_range_core.hpp"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"

#include "molassembler/BondOrders.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/Graph/InnerGraph.h"

#include "temple/Functional.h"

namespace Scine {
namespace molassembler {
namespace interpret {

namespace detail {

Scine::Utils::PositionCollection paste(const std::vector<Scine::Utils::Position>& positions) {
  const unsigned N = positions.size();
  auto matrix = Scine::Utils::PositionCollection(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    matrix.row(i) = positions.at(i);
  }
  return matrix;
}

} // namespace detail

struct MoleculeParts {
  InnerGraph graph;
  std::vector<Scine::Utils::Position> angstromPositions;
  boost::optional<
    std::vector<BondIndex>
  > bondStereopermutatorCandidatesOptional;
};

struct Parts {
  std::vector<MoleculeParts> precursors;
  std::vector<unsigned> componentMap;
  unsigned nZeroLengthPositions = 0;
};

Parts construeParts(
  const Scine::Utils::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const Scine::Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  Parts parts;

  // Discretize bond orders
  const unsigned N = elements.size();

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

  InnerGraph atomCollectionGraph {
    static_cast<InnerGraph::Vertex>(N)
  };

  if(discretization == BondDiscretizationOption::Binary) {
    for(unsigned i = 0; i < N; ++i) {
      for(unsigned j = i + 1; j < N; ++j) {
        double bondOrder = bondOrders.getOrder(i, j);

        if(bondOrder > 0.5) {
          atomCollectionGraph.addEdge(i, j, BondType::Single);
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

          atomCollectionGraph.addEdge(i, j, bond);
        }
      }
    }
  }

  // Calculate the number of components and keep the component map
  unsigned numComponents = atomCollectionGraph.connectedComponents(parts.componentMap);

  parts.precursors.resize(numComponents);

  if(stereopermutatorBondOrderThresholdOptional) {
    for(auto& precursor : parts.precursors) {
      // Empty-vector-initialize the candidate optionals
      precursor.bondStereopermutatorCandidatesOptional = std::vector<BondIndex> {};
    }
  }

  // Map from original index to component index
  std::vector<InnerGraph::Vertex> indexInComponentMap (N);

  /* Maybe
   * - filtered_graph using predicate of componentMap number
   * - copy_graph to new OuterGraph keeping element types and bond orders
   *
   * - alternately, must keep a map of atomcollection index to precursor index
   *   and new precursor atom index in order to transfer edges too
   */

  for(unsigned i = 0; i < N; ++i) {
    auto& precursor = parts.precursors.at(
      parts.componentMap.at(i)
    );

    // Add a new vertex with element information
    InnerGraph::Vertex newIndex = precursor.graph.addVertex(elements.at(i));

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
  for(
    const InnerGraph::Edge& edge :
    boost::make_iterator_range(
      atomCollectionGraph.edges()
    )
  ) {
    InnerGraph::Vertex source = atomCollectionGraph.source(edge),
                       target = atomCollectionGraph.target(edge);

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
      stereopermutatorBondOrderThresholdOptional
      && bondOrders.getOrder(source, target) >= *stereopermutatorBondOrderThresholdOptional
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
  const Scine::Utils::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const Scine::Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  Parts parts = construeParts(
    elements,
    angstromWrapper,
    bondOrders,
    discretization,
    stereopermutatorBondOrderThresholdOptional
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
        OuterGraph {std::move(precursor.graph)},
        AngstromWrapper(detail::paste(precursor.angstromPositions), LengthUnit::Angstrom),
        precursor.bondStereopermutatorCandidatesOptional
      );
    }
  } else {
    result.molecules.reserve(parts.precursors.size());
    for(auto& precursor : parts.precursors) {
      result.molecules.emplace_back(
        OuterGraph {std::move(precursor.graph)}
      );
    }
  }

  result.componentMap = std::move(parts.componentMap);

  return result;
}

MoleculesResult molecules(
  const Scine::Utils::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  return molecules(
    elements,
    angstromWrapper,
    uffBondOrders(elements, angstromWrapper),
    discretization,
    stereopermutatorBondOrderThresholdOptional
  );
}

MoleculesResult molecules(
  const Scine::Utils::AtomCollection& atomCollection,
  const Scine::Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  return molecules(
    atomCollection.getElements(),
    AngstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    discretization,
    stereopermutatorBondOrderThresholdOptional
  );
}

MoleculesResult molecules(
  const Scine::Utils::AtomCollection& atomCollection,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  AngstromWrapper angstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr};

  return molecules(
    atomCollection.getElements(),
    angstromWrapper,
    uffBondOrders(atomCollection.getElements(), angstromWrapper),
    discretization,
    stereopermutatorBondOrderThresholdOptional
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
  const AngstromWrapper& angstromWrapper,
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
    AngstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    discretization
  );
}

} // namespace interpret
} // namespace molassembler
} // namespace Scine
