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

namespace Scine {

namespace molassembler {

namespace detail {

Scine::Utils::PositionCollection paste(const std::vector<Scine::Utils::Position>& positions) {
  const unsigned N = positions.size();
  auto matrix = Scine::Utils::PositionCollection(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    matrix.row(i) = positions.at(i);
  }
  return matrix;
}

} // namespace

InterpretResult interpret(
  const Scine::Utils::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const Scine::Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
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
  std::vector<unsigned> componentMap;
  unsigned numComponents = atomCollectionGraph.connectedComponents(componentMap);

  struct MoleculeParts {
    InnerGraph graph;
    std::vector<Scine::Utils::Position> angstromPositions;
    boost::optional<
      std::vector<BondIndex>
    > bondStereopermutatorCandidatesOptional;
  };

  std::vector<MoleculeParts> moleculePrecursors {numComponents};

  if(stereopermutatorBondOrderThresholdOptional) {
    for(auto& precursor : moleculePrecursors) {
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

  unsigned nZeroLengthPositions = 0;

  for(unsigned i = 0; i < N; ++i) {
    auto& precursor = moleculePrecursors.at(
      componentMap.at(i)
    );

    // Add a new vertex with element information
    InnerGraph::Vertex newIndex = precursor.graph.addVertex(elements.at(i));

    // Save new index in precursor graph
    indexInComponentMap.at(i) = newIndex;

    if(angstromWrapper.positions.row(i).norm() <= 1e-14) {
      ++nZeroLengthPositions;
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
    auto& precursor = moleculePrecursors.at(
      componentMap.at(source)
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

  // Collect results
  InterpretResult result;

  /* Transform precursors into Molecules. Positions may only be used if there
   * is at most one position very close to (0, 0, 0). Otherwise, we assume that
   * the given positions are faulty or no positional information is present,
   * and only the graph is used to create the Molecules.
   */
  if(nZeroLengthPositions < 2) {
    result.molecules.reserve(moleculePrecursors.size());
    for(auto& precursor : moleculePrecursors) {
      result.molecules.emplace_back(
        OuterGraph {std::move(precursor.graph)},
        AngstromWrapper(detail::paste(precursor.angstromPositions), LengthUnit::Angstrom),
        precursor.bondStereopermutatorCandidatesOptional
      );
    }
  } else {
    result.molecules.reserve(moleculePrecursors.size());
    for(auto& precursor : moleculePrecursors) {
      result.molecules.emplace_back(
        OuterGraph {std::move(precursor.graph)}
      );
    }
  }

  result.componentMap = std::move(componentMap);

  return result;
}

InterpretResult interpret(
  const Scine::Utils::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  return interpret(
    elements,
    angstromWrapper,
    uffBondOrders(elements, angstromWrapper),
    discretization,
    stereopermutatorBondOrderThresholdOptional
  );
}

InterpretResult interpret(
  const Scine::Utils::AtomCollection& atomCollection,
  const Scine::Utils::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  return interpret(
    atomCollection.getElements(),
    AngstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    discretization,
    stereopermutatorBondOrderThresholdOptional
  );
}

InterpretResult interpret(
  const Scine::Utils::AtomCollection& atomCollection,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorBondOrderThresholdOptional
) {
  AngstromWrapper angstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr};

  return interpret(
    atomCollection.getElements(),
    angstromWrapper,
    uffBondOrders(atomCollection.getElements(), angstromWrapper),
    discretization,
    stereopermutatorBondOrderThresholdOptional
  );
}

} // namespace molassembler

} // namespace Scine
