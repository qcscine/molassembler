#include "Interpret.h"

#include "BondOrders.h"
#include "boost/graph/connected_components.hpp"

namespace molassembler {

InterpretResult interpret(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const Delib::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization
) {
  // Discretize bond orders
  const int N = elements.size();

  GraphType atomCollectionGraph {
    static_cast<GraphType::vertices_size_type>(N)
  };

  if(discretization == BondDiscretizationOption::Binary) {
    for(int i = 0; i < N; ++i) {
      for(int j = i + 1; j < N; ++j) {
        double bondOrder = bondOrders.getOrder(i, j);

        if(bondOrder > 0.5) {
          auto edgeAddPair = boost::add_edge(i, j, atomCollectionGraph);
          atomCollectionGraph[edgeAddPair.first].bondType = BondType::Single;
        }
      }
    }
  } else if(discretization == BondDiscretizationOption::UFF) {
    for(int i = 0; i < N; ++i) {
      for(int j = i + 1; j < N; ++j) {
        double bondOrder = bondOrders.getOrder(i, j);

        if(bondOrder > 0.5) {
          auto bond = static_cast<BondType>(
            std::round(bondOrder) - 1
          );

          if(bondOrder > 6.5) {
            bond = BondType::Sextuple;
          }

          auto edgeAddPair = boost::add_edge(i, j, atomCollectionGraph);

          atomCollectionGraph[edgeAddPair.first].bondType = bond;
        }
      }
    }
  }

  // Calculate the number of components and keep the component map
  std::vector<unsigned> componentMap(N);

  unsigned numComponents = boost::connected_components(atomCollectionGraph, &componentMap[0]);

  struct MoleculeParts {
    GraphType graph;
    AngstromWrapper angstromWrapper;
  };

  std::vector<MoleculeParts> moleculePrecursors {numComponents};
  std::vector<AtomIndexType> indexInComponentMap (N);

  /* Maybe
   * - filtered_graph using predicate of componentMap number
   * - copy_graph to new GraphType keeping element types and bond orders
   *
   * - alternately, must keep a map of atomcollection index to precursor index
   *   and new precursor atom index in order to transfer edges too
   */

  unsigned nZeroLengthPositions = 0;

  for(int i = 0; i < N; ++i) {
    auto& precursor = moleculePrecursors.at(
      componentMap.at(i)
    );

    AtomIndexType newIndex = boost::add_vertex(precursor.graph);

    // Save new index in precursor graph
    indexInComponentMap.at(i) = newIndex;

    // Copy over the element information into the precursor
    precursor.graph[newIndex].elementType = elements.at(i);

    if(angstromWrapper.positions.at(i).asEigenVector().norm() <= 1e-14) {
      ++nZeroLengthPositions;
    }

    // Copy over position information adjusted by lengthScale
    precursor.angstromWrapper.positions.push_back(
      angstromWrapper.positions.at(i)
    );
  }

  // Copy over edges and bond orders
  for(
    const auto& edge : RangeForTemporary<GraphType::edge_iterator>(
      boost::edges(atomCollectionGraph)
    )
  ) {
    AtomIndexType source = boost::source(edge, atomCollectionGraph),
                  target = boost::target(edge, atomCollectionGraph);

    // Both source and target are part of the same component (since they are bonded)
    auto& precursor = moleculePrecursors.at(
      componentMap.at(source)
    );

    // Copy over the edge
    auto edgeAddPair = boost::add_edge(
      indexInComponentMap.at(source),
      indexInComponentMap.at(target),
      precursor.graph
    );

    precursor.graph[edgeAddPair.first].bondType = atomCollectionGraph[edge].bondType;
  }

  // Collect results
  InterpretResult result;

  /* Transform precursors into Molecules. Positions may only be used if there
   * is at most one position very close to (0, 0, 0). Otherwise, we assume that
   * the given positions are faulty or no positional information is present,
   * and only the graph is used to create the Molecules.
   */
  if(nZeroLengthPositions < 2) {
    result.molecules = temple::map(
      moleculePrecursors,
      [](const MoleculeParts& precursor) -> Molecule {
        return {
          precursor.graph,
          precursor.angstromWrapper
        };
      }
    );
  } else {
    result.molecules = temple::map(
      moleculePrecursors,
      [](const MoleculeParts& precursor) -> Molecule {
        return Molecule {precursor.graph};
      }
    );
  }

  result.componentMap = std::move(componentMap);

  return result;
}

InterpretResult interpret(
  const Delib::ElementTypeCollection& elements,
  const AngstromWrapper& angstromWrapper,
  const BondDiscretizationOption discretization
) {
  return interpret(
    elements,
    angstromWrapper,
    uffBondOrders(elements, angstromWrapper),
    discretization
  );
}

InterpretResult interpret(
  const Delib::AtomCollection& atomCollection,
  const Delib::BondOrderCollection& bondOrders,
  const BondDiscretizationOption discretization
) {
  return interpret(
    atomCollection.getElements(),
    AngstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    discretization
  );
}

InterpretResult interpret(
  const Delib::AtomCollection& atomCollection,
  const BondDiscretizationOption discretization
) {
  AngstromWrapper angstromWrapper {atomCollection.getPositions(), LengthUnit::Bohr};

  return interpret(
    atomCollection.getElements(),
    angstromWrapper,
    uffBondOrders(atomCollection.getElements(), angstromWrapper),
    discretization
  );
}

} // namespace molassembler
