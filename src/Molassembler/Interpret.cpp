/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Interpret.h"

#include "Utils/Geometry/AtomCollection.h"
#include "Utils/Bonds/BondOrderCollection.h"
#include "Utils/Bonds/BondDetector.h"

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
  const AngstromPositions& wrapper,
  const AtomIndex v,
  const std::vector<AtomIndex>& site
) {
  const auto& positions = wrapper.positions;
  const unsigned siteSize = site.size();

  assert(siteSize >= 1 && "Passed hapticPlaneGeometry a non-haptic site!");

  Eigen::Vector3d siteCentroid = Eigen::Vector3d::Zero();
  for(AtomIndex siteAtom : site) {
    siteCentroid += positions.row(siteAtom);
  }
  siteCentroid /= siteSize;

  if(siteSize == 2) {
    const double frontAngle = Cartesian::angle(
      positions.row(v),
      siteCentroid,
      positions.row(site.front())
    );
    const double backAngle = Cartesian::angle(
      positions.row(v),
      siteCentroid,
      positions.row(site.back())
    );

    return HapticPlaneGeometry {
      M_PI / 2.0 - std::min(frontAngle, backAngle),
      0.0
    };
  }

  const Eigen::Vector3d centroidAxis = siteCentroid - positions.row(v).transpose();

  Utils::PositionCollection hapticAtoms(siteSize, 3);
  for(unsigned i = 0; i < siteSize; ++i) {
    hapticAtoms.row(i) = positions.row(site.at(i));
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
  const AngstromPositions& angstromPositions,
  const AtomIndex v
) {
  auto sites = GraphAlgorithms::sites(graph, v);
  if(sites.size() <= 1) {
    return boost::none;
  }

  const auto& positions = angstromPositions.positions;

  const unsigned S = sites.size();
  Eigen::Matrix<double, 3, Eigen::Dynamic> sitePositions(3, S + 1);
  for(unsigned i = 0; i < S; ++i) {
    Eigen::Vector3d averagePosition = Eigen::Vector3d::Zero();
    for(unsigned j : sites[i]) {
      averagePosition += positions.row(j);
    }
    sitePositions.col(i) = averagePosition / sites[i].size();
  }
  sitePositions.col(S) = positions.row(v);

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
  const AngstromPositions& positions,
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
  const unsigned nComponents = countComponents();

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
  const unsigned nComponents = countComponents();
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

  const unsigned N = map.size();
  for(unsigned i = 0; i < N; ++i) {
    const unsigned moleculeIndex = map.at(i);
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

unsigned ComponentMap::countComponents() const {
  assert(!map.empty());
  return *std::max_element(
    std::begin(map),
    std::end(map)
  ) + 1;
}

std::vector<AngstromPositions> ComponentMap::apply(
  const AngstromPositions& positions
) const {
  const unsigned nComponents = countComponents();
  std::vector<unsigned> componentSizes(nComponents, 0);
  for(unsigned i : map) {
    componentSizes.at(i) += 1;
  }

  /* Allocate the component position objects */
  std::vector<AngstromPositions> collections = Temple::map(
    componentSizes,
    [](const unsigned size) { return AngstromPositions(size); }
  );
  std::vector<unsigned> collectionSizeCount(nComponents, 0);

  const unsigned N = map.size();
  for(unsigned i = 0; i < N; ++i) {
    const unsigned moleculeIndex = map.at(i);
    AngstromPositions& componentPositions = collections.at(moleculeIndex);
    unsigned& collectionSize = collectionSizeCount.at(moleculeIndex);
    componentPositions.positions.row(collectionSize) = positions.positions.row(i);
    ++collectionSize;
  }

  return collections;
}

std::vector<PeriodicBoundaryDuplicates> ComponentMap::apply(
  const std::unordered_set<unsigned>& uninterestingAtoms,
  const std::unordered_map<unsigned, unsigned>& ghostAtomMap
) const {
  const unsigned nComponents = countComponents();
  std::vector<PeriodicBoundaryDuplicates> containers(nComponents);

  std::vector<unsigned long> indicesInComponent (nComponents, 0);
  const unsigned N = size();
  for(unsigned i = 0; i < N; ++i) {
    const unsigned component = map.at(i);
    unsigned long& indexInComponent = indicesInComponent.at(component);
    PeriodicBoundaryDuplicates& componentContainers = containers.at(component);

    if(uninterestingAtoms.count(i) > 0) {
      componentContainers.uninterestingAtoms.insert(indexInComponent);
    }

    const auto ghostAtomMapIter = ghostAtomMap.find(i);
    if(ghostAtomMapIter != std::end(ghostAtomMap)) {
      auto transformedMappedIndex = apply(ghostAtomMapIter->second);
      assert(transformedMappedIndex.component == component);
      componentContainers.ghostAtomMap.emplace(
        indexInComponent,
        transformedMappedIndex.atomIndex
      );
    }

    ++indexInComponent;
  }

  return containers;
}

ComponentMap::ComponentIndexPair ComponentMap::apply(const unsigned index) const {
  ComponentIndexPair pair;
  pair.component = map.at(index);
  pair.atomIndex = std::count(
    std::begin(map),
    std::begin(map) + index,
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
  AngstromPositions positions;
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

  /* Must keep a map of atom collection index to precursor index and new
   * precursor atom index in order to transfer edges too
   */

  for(unsigned i = 0; i < N; ++i) {
    auto& precursor = parts.precursors.at(
      parts.componentMap.apply(i).component
    );

    // Add a new vertex with element information
    const PrivateGraph::Vertex newIndex = precursor.graph.addVertex(elements.at(i));

    // Save new index in precursor graph
    indexInComponentMap.at(i) = newIndex;

    if(angstromWrapper.positions.row(i).norm() <= 1e-8) {
      parts.nZeroLengthPositions += 1;
    }
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

  // Split angstrom atom positions and store in precursors
  auto splitPositions = parts.componentMap.apply(angstromWrapper);
  for(unsigned i = 0; i < numComponents; ++i) {
    parts.precursors.at(i).positions = std::move(splitPositions.at(i));
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
        precursor.positions,
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
    Utils::BondDetector::detectBonds(elements, angstromWrapper.getBohr()),
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
    Utils::BondDetector::detectBonds(atomCollection),
    discretization,
    stereopermutatorThreshold
  );
}

/* Periodic variation of molecule instantiation */
MoleculesResult molecules(
  const Utils::AtomCollection& atoms,
  const Utils::BondOrderCollection& periodicBonds,
  const std::unordered_set<unsigned>& uninterestingAtoms,
  const std::unordered_map<unsigned, unsigned>& ghostAtomMap,
  const BondDiscretizationOption discretization,
  const boost::optional<double>& stereopermutatorThreshold
) {
  const AngstromPositions angstroms {atoms.getPositions(), LengthUnit::Bohr};

  Parts parts = construeParts(
    atoms.getElements(),
    angstroms,
    periodicBonds,
    discretization,
    stereopermutatorThreshold
  );

  // Refuse to deal with missing coordinates
  if(parts.nZeroLengthPositions > 1) {
    throw std::runtime_error("Found multiple coordinates of length zero. Please provide atom positions!");
  }

  const auto periodicContainers = parts.componentMap.apply(uninterestingAtoms, ghostAtomMap);

  MoleculesResult result;
  result.molecules.reserve(parts.precursors.size());
  const unsigned N = parts.precursors.size();
  for(unsigned i = 0; i < N; ++i) {
    auto& precursor = parts.precursors.at(i);
    result.molecules.emplace_back(
      Graph {std::move(precursor.graph)},
      precursor.positions,
      precursor.bondStereopermutatorCandidatesOptional,
      periodicContainers.at(i)
    );
  }

  // Copy the component map, removing ghost atoms
  if(ghostAtomMap.empty()) {
    result.componentMap = std::move(parts.componentMap);
  } else {
    for(unsigned i = 0; i < parts.componentMap.size(); ++i) {
      if(ghostAtomMap.count(i) > 0) {
        break;
      }

      result.componentMap.map.push_back(parts.componentMap.map.at(i));
    }
  }

  return result;
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
      auto classificationUncertainty = minimumClassificationProbability(part.graph, part.positions, v);
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

        const auto geometry = hapticPlaneGeometry(part.positions, v, site);

        // Less than 30° and a good plane fit indicate the haptic site is fine
        if(geometry.angle < M_PI / 6 && geometry.rmsd < 0.2) {
          continue;
        }

        if(siteSize == 2) {
          // Suggest the vertex further from the center
          const double frontDistance = (
            part.positions.positions.row(v)
            - part.positions.positions.row(site.front())
          ).norm();
          const double backDistance = (
            part.positions.positions.row(v)
            - part.positions.positions.row(site.back())
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
            part.positions,
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
    Temple::sort(haptics);
    FalsePositive& mostLikely = haptics.back();
    bonds.setOrder(mostLikely.i, mostLikely.j, 0.0);
    haptics = badHapticLigandBonds(atoms, bonds);
  }

  // Then do uncertain bonds
  auto uncertains = uncertainBonds(atoms, bonds);
  while(!uncertains.empty()) {
    Temple::sort(uncertains);
    FalsePositive& mostLikely = uncertains.back();
    bonds.setOrder(mostLikely.i, mostLikely.j, 0.0);
    uncertains = uncertainBonds(atoms, bonds);
  }

  return bonds;
}

} // namespace Interpret
} // namespace Molassembler
} // namespace Scine
