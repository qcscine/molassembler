/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Interpret.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"

#include "Molassembler/BondOrders.h"
#include "Molassembler/Detail/Cartesian.h"
#include "Molassembler/Graph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Shapes/ContinuousMeasures.h"
#include "Molassembler/Temple/Adaptors/AllPairs.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Optionals.h"

#include <Eigen/Geometry>

namespace Scine {
namespace Molassembler {
namespace Interpret {
namespace {

Utils::PositionCollection paste(const std::vector<Utils::Position>& positions) {
  const unsigned N = positions.size();
  auto matrix = Utils::PositionCollection(N, 3);
  for(unsigned i = 0; i < N; ++i) {
    matrix.row(i) = positions.at(i);
  }
  return matrix;
}

double vectorAngle(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
  return std::acos(
    a.dot(b) / (
      a.norm() * b.norm()
    )
  );
}

struct HapticPlaneGeometry {
  double angle;
  double rmsd;
};

HapticPlaneGeometry hapticPlaneGeometry(
  const std::vector<Utils::Position> positions,
  const AtomIndex v,
  const std::vector<AtomIndex>& site
) {
  const unsigned siteSize = site.size();

  assert(siteSize >= 1 && "Passed hapticPlaneGeometry a non-haptic site!");

  Eigen::Vector3d siteCentroid = Eigen::Vector3d::Zero();
  for(AtomIndex siteAtom : site) {
    siteCentroid += positions.at(siteAtom);
  }
  siteCentroid /= siteSize;

  if(siteSize == 2) {
    const double frontAngle = Cartesian::angle(
      positions.at(v),
      siteCentroid,
      positions.at(site.front())
    );
    const double backAngle = Cartesian::angle(
      positions.at(v),
      siteCentroid,
      positions.at(site.back())
    );

    return HapticPlaneGeometry {
      M_PI / 2.0 - std::min(frontAngle, backAngle),
      0.0
    };
  }

  const Eigen::Vector3d centroidAxis = siteCentroid - positions.at(v).transpose();

  Utils::PositionCollection hapticAtoms(siteSize, 3);
  for(unsigned i = 0; i < siteSize; ++i) {
    hapticAtoms.row(i) = positions.at(site.at(i));
  }
  const auto plane = Cartesian::planeOfBestFit(hapticAtoms);
  const double firstAngle = vectorAngle(plane.normal(), centroidAxis);
  const double secondAngle = vectorAngle(plane.normal(), -centroidAxis);
  return HapticPlaneGeometry {
    std::min(firstAngle, secondAngle),
    Cartesian::planeRmsd(plane, hapticAtoms, Temple::iota<AtomIndex>(siteSize))
  };
}

boost::optional<double> minimumClassificationProbability(
  const PrivateGraph& graph,
  const std::vector<Utils::Position>& positions,
  const AtomIndex v
) {
  auto sites = GraphAlgorithms::sites(graph, v);
  if(sites.size() <= 1) {
    return boost::none;
  }

  const unsigned S = sites.size();
  Eigen::Matrix<double, 3, Eigen::Dynamic> sitePositions(3, S + 1);
  for(unsigned i = 0; i < S; ++i) {
    Eigen::Vector3d averagePosition = Eigen::Vector3d::Zero();
    for(unsigned j : sites[i]) {
      averagePosition += positions.at(j);
    }
    sitePositions.col(i) = averagePosition / sites[i].size();
  }
  sitePositions.col(S) = positions.at(v);

  auto normalizedPositions = Shapes::Continuous::normalize(sitePositions);

  // Classify all suitable shapes
  std::vector<Shapes::Shape> viableShapes;
  for(const Shapes::Shape shape : Shapes::allShapes) {
    if(Shapes::size(shape) == S) {
      viableShapes.push_back(shape);
    }
  }

  const auto classifications = Temple::map(viableShapes,
    [&](const Shapes::Shape shape) -> boost::optional<double> {
      const double measure = Shapes::Continuous::shapeCentroidLast(normalizedPositions, shape).measure;
      return Shapes::Continuous::probabilityRandomCloud(measure, shape);
    }
  );

  if(Temple::all_of(classifications)) {
    const double minimumProbability = std::min_element(
      std::begin(classifications),
      std::end(classifications),
      [](const auto& a, const auto& b) -> bool {
        return a.value() < b.value();
      }
    )->value();

    return minimumProbability;
  }

  return boost::none;
}

struct BestRemovals {
  std::vector<AtomIndex> removals;
  double certainty;
};

boost::optional<BestRemovals> bestRemovalFromHapticSite(
  const PrivateGraph& graph,
  const std::vector<Utils::Position>& positions,
  const AtomIndex v,
  const std::vector<std::vector<AtomIndex>>& sites,
  const unsigned hapticSiteIndex
) {
  const unsigned S = sites.size();
  const std::vector<AtomIndex>& hapticSite = sites.at(hapticSiteIndex);
  const unsigned hapticSiteSize = hapticSite.size();
  return Temple::accumulate(
    hapticSite,
    boost::optional<BestRemovals>(boost::none),
    [&](const auto& carry, const AtomIndex siteVertexToRemove) -> boost::optional<BestRemovals> {
      auto graphCopy = graph;
      graphCopy.removeEdge(graphCopy.edge(v, siteVertexToRemove));
      GraphAlgorithms::updateEtaBonds(graphCopy);
      const auto newSites = GraphAlgorithms::sites(graphCopy, v);

      /* Possible effects of bond removal
       * - Separating a haptic ligand into two single-atom ligands,
       *   changing shapes
       *   - recognizable by size change
       *   - accept if new certainty is low or significantly
       *     lower than old (factor 0.5)
       * - Improves the haptic plane angle by removing a bad bond
       *   - recognized by matching sites and comparing angles
       *   - accept if angle halves
       *
       * How to compare cases?
       */
      if(newSites.size() > S) {
        const auto priorCertainty = minimumClassificationProbability(graph, positions, v);
        const auto posteriorCertainty = minimumClassificationProbability(graphCopy, positions, v);

        if(priorCertainty && posteriorCertainty) {
          if(
            posteriorCertainty.value() <= 0.01
            || posteriorCertainty.value() <= 0.5 * priorCertainty.value()
          ) {
            return BestRemovals {{siteVertexToRemove}, 1 - *posteriorCertainty};
          }
        }
      } else {
        const auto siteFindIter = Temple::find_if(
          newSites,
          [&](const auto& newSite) -> bool {
            if(newSite.size() != hapticSiteSize - 1) {
              return false;
            }

            // Match all vertices except the one to remove
            for(const AtomIndex oldSiteAtomIndex : hapticSite) {
              if(oldSiteAtomIndex == siteVertexToRemove) {
                continue;
              }

              if(Temple::find(newSite, oldSiteAtomIndex) == std::end(newSite)) {
                return false;
              }
            }

            return true;
          }
        );

        if(siteFindIter != std::end(newSites)) {
          const double newHapticAngle = hapticPlaneGeometry(positions, v, *siteFindIter).angle;
          const double newHapticAngleCertainty = 1 - newHapticAngle * 2 / M_PI;
          const double carryCertainty = Temple::Optionals::map(carry,
            [](const auto& removal) { return removal.certainty; }
          ).value_or(0.0);
          if(newHapticAngleCertainty > carryCertainty) {
            return BestRemovals {{siteVertexToRemove}, newHapticAngleCertainty};
          }
        }
      }

      return carry;
    }
  );
}

} // namespace

std::vector<std::vector<unsigned>> ComponentMap::invert() const {
  const unsigned nComponents = *std::max_element(
    std::begin(map),
    std::end(map)
  ) + 1;

  std::vector<
    std::vector<unsigned>
  > inverseMaps (nComponents);

  const unsigned N = map.size();
  for(unsigned i = 0; i < N; ++i) {
    inverseMaps.at(map.at(i)).push_back(i);
  }

  return inverseMaps;
}

std::vector<Utils::AtomCollection> ComponentMap::apply(
  const Utils::AtomCollection& atomCollection
) const {
  const unsigned nComponents = *std::max_element(
    std::begin(map),
    std::end(map)
  ) + 1;
  std::vector<unsigned> componentSizes(nComponents, 0);
  for(unsigned i : map) {
    componentSizes.at(i) += 1;
  }

  /* Allocate the collections */
  std::vector<Utils::AtomCollection> collections = Temple::map(
    componentSizes,
    [](const unsigned size) { return Utils::AtomCollection(size); }
  );
  std::vector<unsigned> collectionSizeCount(nComponents, 0);

  for(unsigned i = 0; i < map.size(); ++i) {
    unsigned moleculeIndex = map.at(i);
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

ComponentMap::ComponentIndexPair ComponentMap::apply(const unsigned left) const {
  ComponentIndexPair pair;
  pair.component = map.at(left);
  pair.atomIndex = std::count(
    std::begin(map),
    std::begin(map) + left,
    pair.component
  );
  return pair;
}

unsigned ComponentMap::invert(const ComponentIndexPair& pair) const {
  unsigned count = 0;
  const unsigned N = map.size();
  for(unsigned i = 0; i < N; ++i) {
    if(map[i] == pair.component) {
      ++count;
    }

    if(count == pair.atomIndex + 1) {
      return i;
    }
  }

  throw std::out_of_range("No match found in component map!");
}

struct MoleculeParts {
  PrivateGraph graph;
  std::vector<Utils::Position> angstromPositions;
  boost::optional<
    std::vector<BondIndex>
  > bondStereopermutatorCandidatesOptional;
};

struct Parts {
  std::vector<MoleculeParts> precursors;
  ComponentMap componentMap;
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
  const unsigned numComponents = atomCollectionGraph.connectedComponents(parts.componentMap.map);
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
      parts.componentMap.apply(i).component
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
      parts.componentMap.apply(source).component
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

std::vector<FalsePositive> uncertainBonds(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders
) {
  std::vector<FalsePositive> falsePositives;

  Parts parts = construeParts(
    atomCollection.getElements(),
    AngstromPositions {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    BondDiscretizationOption::Binary,
    boost::none
  );

  for(unsigned component = 0; component < parts.precursors.size(); ++component) {
    MoleculeParts& part = parts.precursors[component];
    GraphAlgorithms::updateEtaBonds(part.graph);

    std::vector<std::pair<PrivateGraph::Vertex, double>> sketchyClassifications;

    for(const PrivateGraph::Vertex v : part.graph.vertices()) {
      /* Detect high uncertainty shape classifications. Adjacent pairs with high
       * uncertainties are candidates for false positives.
       */
      auto classificationUncertainty = minimumClassificationProbability(part.graph, part.angstromPositions, v);
      if(classificationUncertainty && *classificationUncertainty >= 0.5) {
        sketchyClassifications.emplace_back(v, *classificationUncertainty);
      }
    }

    // NOTE: Cannot remove bonds with overlapping constituent atoms!
    for(const auto& sketchyPair : Temple::Adaptors::allPairs(sketchyClassifications)) {
      const unsigned i = sketchyPair.first.first;
      const unsigned j = sketchyPair.second.first;
      const double i_random = sketchyPair.first.second;
      const double j_random = sketchyPair.second.second;

      if(part.graph.edgeOption(i, j)) {
        falsePositives.push_back(
          FalsePositive {
            parts.componentMap.invert(ComponentMap::ComponentIndexPair {component, i}),
            parts.componentMap.invert(ComponentMap::ComponentIndexPair {component, j}),
            i_random * j_random
          }
        );
      }
    }
  }

  return falsePositives;
}

std::vector<FalsePositive> badHapticLigandBonds(
  const Utils::AtomCollection& atomCollection,
  const Utils::BondOrderCollection& bondOrders
) {
  std::vector<FalsePositive> falsePositives;

  // Avoid duplicate bonds in list with reverse i-j
  auto addFalsePositive = [&](unsigned i, unsigned j, double p) {
    if(i > j) {
      std::swap(i, j);
    }
    auto findIter = std::find_if(
      std::begin(falsePositives),
      std::end(falsePositives),
      [=](const FalsePositive& fp) -> bool {
        return fp.i == i && fp.j == j;
      }
    );
    if(findIter == std::end(falsePositives)) {
      falsePositives.push_back(FalsePositive {i, j, p});
    }
  };

  Parts parts = construeParts(
    atomCollection.getElements(),
    AngstromPositions {atomCollection.getPositions(), LengthUnit::Bohr},
    bondOrders,
    BondDiscretizationOption::Binary,
    boost::none
  );

  for(unsigned component = 0; component < parts.precursors.size(); ++component) {
    MoleculeParts& part = parts.precursors[component];
    GraphAlgorithms::updateEtaBonds(part.graph);

    for(const PrivateGraph::Vertex v : part.graph.vertices()) {
      /* Detect haptic shape planes with large angles to the axis defined by the
       * site position and the central atom or with high rms deviations on
       * their plane fit
       */
      const auto sites = GraphAlgorithms::sites(part.graph, v);
      const unsigned S = sites.size();
      for(unsigned siteIndex = 0; siteIndex < S; ++siteIndex) {
        const auto& site = sites.at(siteIndex);
        const unsigned siteSize = site.size();
        if(siteSize == 1) {
          continue;
        }

        const auto geometry = hapticPlaneGeometry(part.angstromPositions, v, site);

        // Less than 30° and a good plane fit indicate the haptic site is fine
        if(geometry.angle < M_PI / 6 && geometry.rmsd < 0.2) {
          continue;
        }

        if(siteSize == 2) {
          // Suggest the vertex further from the center
          const double frontDistance = (
            part.angstromPositions.at(v)
            - part.angstromPositions.at(site.front())
          ).norm();
          const double backDistance = (
            part.angstromPositions.at(v)
            - part.angstromPositions.at(site.back())
          ).norm();
          const AtomIndex toRemove = frontDistance < backDistance ? site.back() : site.front();
          addFalsePositive(
            parts.componentMap.invert(
              ComponentMap::ComponentIndexPair {component, v}
            ),
            parts.componentMap.invert(
              ComponentMap::ComponentIndexPair {component, toRemove}
            ),
            geometry.angle * 2 / M_PI
          );
        } else {
          const auto suggestedRemovalOption = bestRemovalFromHapticSite(
            part.graph,
            part.angstromPositions,
            v,
            sites,
            siteIndex
          );
          if(suggestedRemovalOption) {
            for(const AtomIndex w : suggestedRemovalOption->removals) {
              addFalsePositive(
                parts.componentMap.invert(
                  ComponentMap::ComponentIndexPair {component, v}
                ),
                parts.componentMap.invert(
                  ComponentMap::ComponentIndexPair {component, w}
                ),
                suggestedRemovalOption->certainty
              );
            }
          }
        }
      }
    }
  }
  return falsePositives;
}

Utils::BondOrderCollection removeFalsePositives(
  const Utils::AtomCollection& atoms,
  Utils::BondOrderCollection bonds
) {
  // First do bad haptic bond orders
  auto haptics = badHapticLigandBonds(atoms, bonds);
  while(!haptics.empty()) {
    Temple::InPlace::sort(haptics);
    FalsePositive& mostLikely = haptics.back();
    bonds.setOrder(mostLikely.i, mostLikely.j, 0.0);
    haptics = badHapticLigandBonds(atoms, bonds);
  }

  // Then do uncertain bonds
  auto uncertains = uncertainBonds(atoms, bonds);
  while(!uncertains.empty()) {
    Temple::InPlace::sort(uncertains);
    FalsePositive& mostLikely = uncertains.back();
    bonds.setOrder(mostLikely.i, mostLikely.j, 0.0);
    uncertains = uncertainBonds(atoms, bonds);
  }

  return bonds;
}

} // namespace Interpret
} // namespace Molassembler
} // namespace Scine
