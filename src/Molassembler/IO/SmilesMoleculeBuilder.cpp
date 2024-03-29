/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "Molassembler/IO/SmilesMoleculeBuilder.h"

#include "Molassembler/IO/SmilesCommon.h"
#include "Molassembler/IO/SmilesBondStereo.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/Stereopermutators/AbstractPermutations.h"
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"
#include "Molassembler/Stereopermutators/ShapeVertexMaps.h"
#include "Molassembler/RankingInformation.h"
#include "Molassembler/Graph.h"
#include "Molassembler/Graph/PrivateGraph.h"
#include "Molassembler/Graph/GraphAlgorithms.h"
#include "Molassembler/Modeling/BondDistance.h"
#include "Molassembler/Modeling/ShapeInference.h"
#include "Molassembler/Molecule.h"
#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Shapes/Properties.h"

#include "Molassembler/Temple/Adaptors/Combinations.h"
#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Temple/Functional.h"

#include "boost/graph/max_cardinality_matching.hpp"

#include <iostream>
#include <fstream>

namespace Scine {
namespace Molassembler {
namespace IO {

bool PiSubgraph::hasUnpairedElectrons(
  const Vertex i,
  const int charge,
  const PrivateGraph& g
) {
  const Utils::ElementType e = g.elementType(i);
  const int electrons = Utils::ElementInfo::valElectrons(e);
  const int valence = vertexValence(i, g);
  return (electrons - valence - charge) > 0;
}

boost::optional<PiSubgraph::EligibleOmissible> PiSubgraph::multipleOrderAdjacent(
  Vertex i,
  const PrivateGraph& component,
  const AtomData& atomData,
  const boost::optional<unsigned>& neighborCount
) {
  const Utils::ElementType e = Utils::ElementInfo::base(component.elementType(i));

  const bool isMultipleOrderAdjacent = Temple::any_of(
    component.edges(i),
    [&](const PrivateGraph::Edge edge) -> bool {
      const BondType b = component.bondType(edge);
      return b == BondType::Double || b == BondType::Triple;
    }
  );

  if(!isMultipleOrderAdjacent) {
    return boost::none;
  }

  // Exceptions!
  const bool threeNeighborPositiveNitrogen = (
    e == Utils::ElementType::N
    && atomData.chargeOptional == 1
    && neighborCount == 3U
  );
  const bool unchargedSulfur = (
    e == Utils::ElementType::S
    && !atomData.chargeOptional
  );

  if(!threeNeighborPositiveNitrogen && !unchargedSulfur) {
    // Neither exception applies, so non-eligible and non-omissible
    return EligibleOmissible {false, false};
  }

  if(unchargedSulfur) {
    /* Sulfur is valence-filled, so as long as the neighbor count isn't
     * fixed or is exactly four, it's eligible. It's omissible only if the
     * neighbor count is variable.
     */
    return EligibleOmissible {
      !neighborCount || neighborCount == 4U,
      !neighborCount
    };
  }

  // Three-neighbor positive nitrogen is the remaining case
  return EligibleOmissible {true, false};
}

boost::optional<PiSubgraph::EligibleOmissible> PiSubgraph::threeNeighborChargedCarbon(
  Vertex i,
  const PrivateGraph& component,
  const AtomData& atomData,
  const boost::optional<unsigned>& neighborCount
) {
  if(
    Utils::ElementInfo::base(component.elementType(i)) == Utils::ElementType::C
    && Temple::Optionals::map(
      atomData.chargeOptional,
      [](int x) { return std::abs(x); }
    ) == 1
    && neighborCount == 3U
  ) {
    return EligibleOmissible {false, false};
  }

  return boost::none;
}

boost::optional<PiSubgraph::EligibleOmissible> PiSubgraph::neutralTrivalents(
  Vertex i,
  const PrivateGraph& component,
  const AtomData& atomData,
  const boost::optional<unsigned>& neighborCount
) {
  const Utils::ElementType e = Utils::ElementInfo::base(component.elementType(i));
  const bool rightElementType = (
    e == Utils::ElementType::N
    || e == Utils::ElementType::P
    || e == Utils::ElementType::As
    || e == Utils::ElementType::Sb
  );

  if(!rightElementType) {
    return boost::none;
  }

  const bool twoNeighborNegativeTrivalents = (
    atomData.chargeOptional == -1
    && neighborCount == 2U
  );
  if(twoNeighborNegativeTrivalents) {
    return EligibleOmissible {false, false};
  }
  if(!atomData.chargeOptional) {
    // Eligible if three neighbors or fillable to three neighbors
    return EligibleOmissible {neighborCount < 3U, !neighborCount};
  }

  return boost::none;
}

boost::optional<PiSubgraph::EligibleOmissible> PiSubgraph::neutralDivalents(
  Vertex i,
  const PrivateGraph& component,
  const AtomData& /* atomData */,
  const boost::optional<unsigned>& neighborCount
) {
  const Utils::ElementType e = Utils::ElementInfo::base(component.elementType(i));

  const bool rightElementType = (
    e == Utils::ElementType::O
    || e == Utils::ElementType::S
    || e == Utils::ElementType::Se
    || e == Utils::ElementType::Te
  );

  if(rightElementType) {
    return EligibleOmissible {neighborCount < 2U, !neighborCount};
  }

  return boost::none;
}

PiSubgraph::EligibleOmissible PiSubgraph::eligibleOmissible(
  const Vertex i,
  const PrivateGraph& component,
  const AtomData& atomData
) {
  const Utils::ElementType e = Utils::ElementInfo::base(component.elementType(i));

  if(!PiSubgraph::permittedElementType(e)) {
    return {false, false};
  }

  /* Do not add:
   * - double or triple bond-adjacent atoms
   *   - unless charged nitrogen with three neighbors or uncharged sulphur with four
   *     (note sulphur can be valence filled and is hence omissible if not a bracket atom)
   * - carbon ions with three neighbors
   *   - ions must be indicated with bracket atoms and are hence not valence filled
   * - N, P, As, Sb
   *   - with three neighbors and neutral
   *   - or with two neighbors and negative charge (bracket atom)
   * - O, S, Se, Te with two neighbors
   *
   * Optionally omissible:
   * - Every eligible atom where valence filling could reach the neighbors requirement
   */
  const auto neighborCount = [&]() -> boost::optional<unsigned> {
    if(atomData.atomBracket || !isValenceFillElement(e)) {
      return component.degree(i) + atomData.hCount.value_or(0);
    }

    /* Valence-fill elements with full valences have fully determined neighbor
     * counts, too. No valence-fill elements have valid valencies in single
     * increments. E.g. for nitrogen 3 and 5. Aromaticity would increment the
     * valence only by one.
     */
    if(isValenceFillElement(e)) {
      const int valence = vertexValence(i, component);
      if(valenceFillElementImplicitHydrogenCount(valence, e) == 0) {
        return component.degree(i);
      }
    }

    return boost::none;
  }();

  if(auto opt = multipleOrderAdjacent(i, component, atomData, neighborCount)) {
    return opt.value();
  }

  if(auto opt = threeNeighborChargedCarbon(i, component, atomData, neighborCount)) {
    return opt.value();
  }

  if(auto opt = neutralTrivalents(i, component, atomData, neighborCount)) {
    return opt.value();
  }

  if(auto opt = neutralDivalents(i, component, atomData, neighborCount)) {
    return opt.value();
  }

  return {true, false};
}

boost::optional<PiSubgraph::VertexSet> PiSubgraph::match() {
  const auto V = boost::num_vertices(graph);

  auto tryMatch = [](const Graph& g) -> bool {
    const unsigned S = boost::num_vertices(g);
    std::vector<PrivateGraph::Vertex> mate(S);
    const bool success = boost::checked_edmonds_maximum_cardinality_matching(g, &mate[0]);
    if(!success) {
      return false;
    }
    const unsigned matchingSize = boost::matching_size(g, &mate[0]);
    return 2 * matchingSize == S;
  };

  // Try to match the full pi subgraph without omissions
  if(tryMatch(graph)) {
    // All vertices are matched, no omissions
    std::unordered_set<Vertex> outerVertices;
    for(auto mapPair : index.left) {
      outerVertices.insert(mapPair.first);
    }
    return outerVertices;
  }

  // Minimum size the subgraph has be is two
  const unsigned maximumOmissions = std::max(V - 2, omissible.size());
  for(unsigned omissions = 1; omissions <= maximumOmissions; ++omissions) {
    for(const auto omit : Temple::Adaptors::combinations(omissible, omissions)) {
      /* TODO This is liable to explode combinatorially and you could use
       * additional criteria to figure out if this choice of omissions is
       * likely to improve anything using e.g. connected components of the pi
       * subgraph.
       */

      auto& subgraph = graph.create_subgraph();
      for(unsigned i = 0; i < V; ++i) {
        if(Temple::find(omit, i) == std::end(omit)) {
          boost::add_vertex(i, subgraph);
        }
      }

      if(tryMatch(subgraph)) {
        std::unordered_set<Vertex> outerVertices;
        const unsigned S = boost::num_vertices(subgraph);
        for(unsigned i = 0; i < S; ++i) {
          const unsigned graphVertex = subgraph.local_to_global(i);
          outerVertices.insert(index.right.at(graphVertex));
        }
        return outerVertices;
      }
    }
  }

  return boost::none;
}

BondType MoleculeBuilder::mutualBondType(
  const boost::optional<BondType>& a,
  const boost::optional<BondType>& b
) {
  /* Ensure that the specified bond order matches. Both bond orders
   * are optionals. If one of both is specified, that one's type is used.
   * If neither is specified, use a single bond. If both are specified,
   * their type is used if it matches, otherwise throw.
   */
  if(!a && !b) {
    return BondType::Single;
  }

  if(a && !b) {
    return a.value();
  }

  if(!a && b) {
    return b.value();
  }

  if(a.value() != b.value()) {
    throw std::runtime_error("Mismatched ring closing bond order");
  }

  return a.value();
}

inline Shapes::Vertex operator "" _v(unsigned long long v) {
  return Shapes::Vertex(v);
}

std::vector<Shapes::Vertex> MoleculeBuilder::shapeMap(const ChiralData& chiralData) {
  switch(chiralData.shape) {
    case Shapes::Shape::Tetrahedron: switch(chiralData.chiralIndex) {
      case 1: return {{0_v, 1_v, 2_v, 3_v}}; // @, TH1
      case 2: return {{0_v, 1_v, 3_v, 2_v}}; // @@, TH2
      default: throw std::out_of_range("No such stereo index for tetrahedron");
    }
    case Shapes::Shape::Square: switch(chiralData.chiralIndex) {
      case 1: return {{0_v, 1_v, 2_v, 3_v}}; // SP1 = U
      case 2: return {{0_v, 2_v, 3_v, 1_v}}; // SP2 = 4
      case 3: return {{3_v, 2_v, 0_v, 1_v}}; // SP3 = Z
      default: throw std::out_of_range("No such stereo index for square");
    }
    case Shapes::Shape::TrigonalBipyramid: switch(chiralData.chiralIndex) {
      case  1: return {{1_v, 2_v, 3_v, 0_v, 4_v}}; // TB1 = a, e, @
      case  2: return {{1_v, 3_v, 2_v, 0_v, 4_v}}; // TB2 = a, e, @@
      case  3: return {{1_v, 2_v, 4_v, 0_v, 3_v}}; // TB3 = a, d, @
      case  4: return {{1_v, 4_v, 2_v, 0_v, 3_v}}; // TB4 = a, d, @@
      case  5: return {{1_v, 3_v, 4_v, 0_v, 2_v}}; // TB5 = a, c, @
      case  6: return {{1_v, 4_v, 3_v, 0_v, 2_v}}; // TB6 = a, c, @@
      case  7: return {{2_v, 3_v, 4_v, 0_v, 1_v}}; // TB7 = a, b, @
      case  8: return {{2_v, 4_v, 3_v, 0_v, 1_v}}; // TB8 = a, b, @@
      case  9: return {{0_v, 2_v, 3_v, 1_v, 4_v}}; // TB9 = b, e, @
      case 10: return {{0_v, 2_v, 4_v, 1_v, 3_v}}; // TB10 = b, d, @
      case 11: return {{0_v, 3_v, 2_v, 1_v, 4_v}}; // TB11 = b, e, @@
      case 12: return {{0_v, 4_v, 2_v, 1_v, 3_v}}; // TB12 = b, d, @@
      case 13: return {{0_v, 3_v, 4_v, 1_v, 2_v}}; // TB13 = b, c, @
      case 14: return {{0_v, 4_v, 3_v, 1_v, 2_v}}; // TB14 = b, c, @@
      case 15: return {{0_v, 1_v, 3_v, 2_v, 4_v}}; // TB15 = c, e, @
      case 16: return {{0_v, 1_v, 4_v, 2_v, 3_v}}; // TB16 = c, d, @
      case 17: return {{0_v, 1_v, 2_v, 3_v, 4_v}}; // TB17 = d, e, @
      case 18: return {{0_v, 2_v, 1_v, 3_v, 4_v}}; // TB18 = d, e, @@
      case 19: return {{0_v, 4_v, 1_v, 2_v, 3_v}}; // TB19 = c, d, @@
      case 20: return {{0_v, 3_v, 1_v, 2_v, 4_v}}; // TB20 = c, e, @@

      default: throw std::out_of_range("No such stereo index for trigonal bipyramid");
    }
    case Shapes::Shape::Octahedron: switch(chiralData.chiralIndex) {
      /* Look along an axis, what remains is a square with a winding. So
       * square shapes are reused with definitions of square shapes and windings
       */
      case  1: return {{1_v, 2_v, 3_v, 4_v, 0_v, 5_v}}; // OH1 = a, f, U, @
      case  2: return {{4_v, 3_v, 2_v, 1_v, 0_v, 5_v}}; // OH2 = a, f, U, @@
      case  3: return {{1_v, 2_v, 3_v, 5_v, 0_v, 4_v}}; // OH3 = a, e, U, @
      case 16: return {{5_v, 3_v, 2_v, 1_v, 0_v, 4_v}}; // OH16 = a, e, U, @@
      case  6: return {{1_v, 2_v, 4_v, 5_v, 0_v, 3_v}}; // OH6 = a, d, U, @
      case 18: return {{5_v, 4_v, 2_v, 1_v, 0_v, 3_v}}; // OH18 = a, d, U, @@
      case 19: return {{1_v, 3_v, 4_v, 5_v, 0_v, 2_v}}; // OH19 = a, c, U, @
      case 24: return {{5_v, 4_v, 3_v, 1_v, 0_v, 2_v}}; // OH24 = a, c, U, @@
      case 25: return {{2_v, 3_v, 4_v, 5_v, 0_v, 1_v}}; // OH25 = a, b, U, @
      case 30: return {{5_v, 4_v, 3_v, 2_v, 0_v, 1_v}}; // OH30 = a, b, U, @@

      /* Note: For the Z shape, the connection between the first two atoms
       * determines the winding.
       */
      case  4: return {{1_v, 2_v, 4_v, 3_v, 0_v, 5_v}}; // OH4 = a, f, Z, @
      case 14: return {{3_v, 4_v, 2_v, 1_v, 0_v, 5_v}}; // OH14 = a, f, Z, @@
      case  5: return {{1_v, 2_v, 5_v, 3_v, 0_v, 4_v}}; // OH5 = a, e, Z, @
      case 15: return {{3_v, 5_v, 2_v, 1_v, 0_v, 4_v}}; // OH15 = a, e, Z, @@
      case  7: return {{1_v, 2_v, 5_v, 4_v, 0_v, 3_v}}; // OH7 = a, d, Z, @
      case 17: return {{4_v, 5_v, 2_v, 1_v, 0_v, 3_v}}; // OH17 = a, d, Z, @@
      case 20: return {{1_v, 3_v, 5_v, 4_v, 0_v, 2_v}}; // OH20 = a, c, Z, @
      case 23: return {{4_v, 5_v, 3_v, 1_v, 0_v, 2_v}}; // OH23 = a, c, Z, @@
      case 26: return {{2_v, 3_v, 5_v, 4_v, 0_v, 1_v}}; // OH26 = a, b, Z, @
      case 29: return {{4_v, 5_v, 3_v, 2_v, 0_v, 1_v}}; // OH29 = a, b, Z, @@

      /* For the 4 shape, the connection between the second and third atom
       * determines the winding.
       */
      case 10: return {{4_v, 2_v, 3_v, 1_v, 0_v, 5_v}}; // OH10 = a, f, 4, @
      case  8: return {{1_v, 3_v, 2_v, 4_v, 0_v, 5_v}}; // OH8 = a, f, 4, @@
      case 11: return {{5_v, 2_v, 3_v, 1_v, 0_v, 4_v}}; // OH11 = a, e, 4, @
      case  9: return {{1_v, 3_v, 2_v, 5_v, 0_v, 4_v}}; // OH9 = a, e, 4, @@
      case 13: return {{5_v, 2_v, 4_v, 1_v, 0_v, 3_v}}; // OH13 = a, d, 4, @
      case 12: return {{1_v, 4_v, 2_v, 5_v, 0_v, 3_v}}; // OH12 = a, d, 4, @@
      case 22: return {{5_v, 3_v, 4_v, 1_v, 0_v, 2_v}}; // OH22 = a, c, 4, @
      case 21: return {{1_v, 4_v, 3_v, 5_v, 0_v, 2_v}}; // OH21 = a, c, 4, @@
      case 28: return {{5_v, 3_v, 4_v, 2_v, 0_v, 1_v}}; // OH28 = a, b, 4, @
      case 27: return {{2_v, 4_v, 3_v, 5_v, 0_v, 1_v}}; // OH27 = a, b, 4, @@

      default: throw std::out_of_range("No such stereo index for octahedron");
    }
    default: throw std::out_of_range("No shape map for selected shape");
  }
}

void MoleculeBuilder::addAtom(AtomData atom) {
  const PrivateGraph::Vertex newVertex = graph.addVertex(atom.getElement());

  if(atom.partialElement.Z == 1 && atom.hCount && atom.hCount.value() != 0) {
    throw std::runtime_error("Hydrogen atoms cannot have hydrogen counts!");
  }

  if(lastBondData) {
    assert(!vertexStack.empty());
    const PrivateGraph::Vertex parent = vertexStack.top();
    const BondData& bond = lastBondData.value();

    // Add edge to graph
    const BondType type = Temple::Optionals::map(bond.type, BondData::toBondType).value_or(BondType::Single);
    graph.addEdge(parent, newVertex, type);

    const bool inferredAromaticBond = (
      !bond.type
      && atom.partialElement.aromatic
      && vertexData.at(parent).partialElement.aromatic
    );

    if(inferredAromaticBond || bond.type == SmilesBondType::Aromatic) {
      piSubgraphEdges.push_back(graph.edge(parent, newVertex));
    }

    // Store stereo-marked bonds for later
    if(bond.type == SmilesBondType::Forward || bond.type == SmilesBondType::Backward) {
      stereoMarkedBonds.emplace_back(
        parent,
        newVertex,
        bond.type.value()
      );
    }
  }

  if(vertexStack.empty()) {
    vertexStack.push(newVertex);
  } else {
    vertexStack.top() = newVertex;
  }

  vertexData.push_back(std::move(atom));
  lastBondData = BondData {};
}

void MoleculeBuilder::addRingClosure(const BondData& bond) {
  assert(bond.ringNumber);
  const auto maybeBondType = Temple::Optionals::map(bond.type, BondData::toBondType);
  const unsigned key = bond.ringNumber.value();
  const auto findIter = ringClosures.find(key);
  if(findIter == std::end(ringClosures)) {
    // Add the entry to the map for later
    ringClosures.emplace_hint(
      findIter,
      std::piecewise_construct,
      std::make_tuple(key),
      std::make_tuple(vertexStack.top(), maybeBondType)
    );
  } else {
    // Add the edge to the graph now and remove the map entry
    const PrivateGraph::Vertex a = findIter->second.first;
    const PrivateGraph::Vertex b = vertexStack.top();

    if(a == b) {
      throw std::runtime_error("Same-atom ring-closing bond!");
    }

    if(graph.edgeOption(a, b) != boost::none) {
      throw std::runtime_error("Ring closing bond already exists!");
    }

    // Ensure the specified bond types match (this fn can throw)
    const BondType mutualType = mutualBondType(
      findIter->second.second,
      maybeBondType
    );

    graph.addEdge(a, b, mutualType);

    const bool inferredAromaticBond = (
      !findIter->second.second
      && !bond.type
      && vertexData.at(a).partialElement.aromatic
      && vertexData.at(b).partialElement.aromatic
    );
    if(inferredAromaticBond || bond.type == SmilesBondType::Aromatic) {
      piSubgraphEdges.push_back(graph.edge(a, b));
    }

    // Remove the entry from the map, freeing up the index
    ringClosures.erase(findIter);
  }
}

void MoleculeBuilder::setShapes(
  std::vector<Molecule>& molecules,
  const std::vector<unsigned>& componentMap,
  const std::vector<PrivateGraph::Vertex>& indexInComponentMap
) {
  const unsigned N = vertexData.size();
  for(unsigned i = 0; i < N; ++i) {
    AtomData& atomData = vertexData.at(i);
    Molecule& mol = molecules.at(componentMap.at(i));
    const AtomIndex atomIndex = indexInComponentMap.at(i);

    if(atomData.chiralOptional) {
      ChiralData& chiralData = atomData.chiralOptional.value();
      const auto& stereoOption = mol.stereopermutators().option(atomIndex);
      if(!stereoOption) {
        continue;
      }
      const unsigned numSites = stereoOption->getRanking().sites.size();
      if(numSites == Shapes::size(chiralData.shape)) {
        mol.setShapeAtAtom(atomIndex, chiralData.shape);
      } else if(chiralData.chiralIndex <= 2) {
        /* We interpreted an @/@@ primarily as Tetrahedron 1 / 2, but
         * they're also in use for trigonal biypramid and octahedron markers
         * as shortcuts for TB1/TB2 and OH1/OH2.
         */
        if(numSites == 5) {
          chiralData.shape = Shapes::Shape::TrigonalBipyramid;
          mol.setShapeAtAtom(atomIndex, chiralData.shape);
        } else if(numSites == 6) {
          chiralData.shape = Shapes::Shape::Octahedron;
          mol.setShapeAtAtom(atomIndex, chiralData.shape);
        }
      }
    } else if(atomData.chargeOptional) {
      /* We need to call VSEPR with the supplied formal charge and make sure
       * the right shape has been inferred. If not, we need to set it.
       */
      const auto& stereoOption = mol.stereopermutators().option(atomIndex);
      if(!stereoOption) {
        // This can happen for ions, for instance, so we consider it harmless
        continue;
      }

      const auto modelArgs = ShapeInference::reduceToSiteInformation(
        mol.graph(),
        atomIndex,
        stereoOption->getRanking()
      );

      const auto shapeOption = ShapeInference::vsepr(
        mol.graph().elementType(atomIndex),
        modelArgs,
        *atomData.chargeOptional
      );

      if(shapeOption && shapeOption.value() != stereoOption->getShape()) {
        mol.setShapeAtAtom(atomIndex, *shapeOption);
      }
    }
  }
}

void MoleculeBuilder::setAtomStereo(
  std::vector<Molecule>& molecules,
  const std::vector<unsigned>& componentMap,
  const std::vector<PrivateGraph::Vertex>& indexInComponentMap,
  const std::string& smiles
) {
  // Aromatic shapes are bent and equilateral triangle only
  const std::map<unsigned, Shapes::Shape> aromaticShapes {
    {2, Shapes::Shape::Bent},
    {3, Shapes::Shape::EquilateralTriangle},
  };

  const unsigned N = vertexData.size();
  for(unsigned i = 0; i < N; ++i) {
    const AtomData& atomData = vertexData.at(i);

    // Set aromatic shapes
    if(atomData.partialElement.aromatic) {
      Molecule& mol = molecules.at(componentMap.at(i));
      if(auto permutator = mol.stereopermutators().option(indexInComponentMap.at(i))) {
        const unsigned size = permutator->getRanking().sites.size();
        if(aromaticShapes.count(size) == 0) {
          throw std::logic_error("Atom " + std::to_string(permutator->placement()) + " with " + std::to_string(size) + " substitutents cannot be aromatic");
        }
        const Shapes::Shape expectedShape = aromaticShapes.at(size);
        if(expectedShape != permutator->getShape()) {
          mol.setShapeAtAtom(indexInComponentMap.at(i), expectedShape);
        }
      }
    }

    if(!atomData.chiralOptional) {
      continue;
    }
    const ChiralData& chiralData = *atomData.chiralOptional;

    /* Set the shape given in the chiral data */
    Molecule& mol = molecules.at(componentMap.at(i));
    auto stereopermutatorOptional = mol.stereopermutators().option(
      indexInComponentMap.at(i)
    );
    if(!stereopermutatorOptional) {
      throw std::logic_error("Atom stereopermutator missing for stereomarked atom!");
    }
    if(stereopermutatorOptional->getShape() != chiralData.shape) {
      throw std::logic_error("Mismatched shape for set chiral data");
    }
    const AtomStereopermutator& permutator = *stereopermutatorOptional;
    if(permutator.numAssignments() < 2) {
      std::cerr << "Warning: Smiles '" << smiles << "' contains a stereo marker for a non-stereogenic " << Shapes::name(chiralData.shape) << " shape center on atom " << permutator.placement() << "\n";
      continue;
    }

    /* Now for the shape stereo markers:
     * - Ordering the sites of the ranking by their constituting indices
     *   yields the order in which they were specified in the SMILES string
     * - We apply any weirdness (like that hcounts get placed at the front of
     *   the list or TODO something about ring closing bonds)
     * - We transfer them onto shape vertex indices depending on the shape and
     *   specified chiral index using shapeMap
     * - We generate a stereopermutation from the siteToShapeVertexMap
     * - And then go looking for it in the list of feasibles
     */
    const unsigned S = Shapes::size(chiralData.shape);
    const RankingInformation& ranking = permutator.getRanking();
    std::vector<SiteIndex> sortedSites = Temple::sorted(
      Temple::iota<SiteIndex>(S),
      [&](const SiteIndex a, const SiteIndex b) -> bool {
        return ranking.sites.at(a) < ranking.sites.at(b);
      }
    );

    /* Atom bracket hcount special case:
     *
     * If one of the neighbors is a hydrogen atom and is represented as an
     * hcount instead of explicitly, then it is considered to be the first
     * atom in the clockwise or anticlockwise counting.
     */
    if(vertexData.at(i).hCount == 1U) {
      for(SiteIndex j {0}; j < ranking.sites.size(); ++j) {
        if(
          ranking.sites.at(j).size() == 1
          && mol.graph().elementType(ranking.sites.at(j).front()) == Utils::ElementType::H
        ) {
          auto explicitHydrogenIter = std::find(std::begin(sortedSites), std::end(sortedSites), j);
          if(explicitHydrogenIter == std::end(sortedSites)) {
            throw std::runtime_error("Failed to find explicit hydrogen site in sorted sites");
          }
          /* Rotate the hcount site index to the first position in the sorted
           * sites
           */
          std::rotate(
            std::begin(sortedSites),
            explicitHydrogenIter,
            explicitHydrogenIter + 1
          );
          break;
        }
      }
    }

    // TODO missing weirdness re: ring closing bonds

    /* Transfer the sorted sites onto shape vertices */
    auto vertexMap = Shapes::Properties::inverseRotation(shapeMap(chiralData));

    std::vector<Shapes::Vertex> siteToShapeVertex(S);
    for(unsigned j = 0; j < S; ++j) {
      siteToShapeVertex.at(sortedSites.at(j)) = vertexMap.at(j);
    }
    auto siteToShapeVertexMap = SiteToShapeVertexMap::from(siteToShapeVertex);

    /* Create a stereopermutation and look for it in the list of feasibles */
    auto soughtStereopermutation = stereopermutationFromSiteToShapeVertexMap(
      siteToShapeVertexMap,
      ranking.links,
      permutator.getAbstract().canonicalSites
    );

    auto soughtRotations = Stereopermutations::generateAllRotations(soughtStereopermutation, chiralData.shape);
    const auto& assignables = permutator.getFeasible();
    auto assignmentIter = Temple::find_if(
      assignables,
      [&](const unsigned stereopermutationIndex) -> bool {
        const auto& stereopermutation = permutator.getAbstract().permutations.list.at(stereopermutationIndex);
        return Temple::find(soughtRotations, stereopermutation) != std::end(soughtRotations);
      }
    );

    if(assignmentIter == std::end(assignables)) {
      throw std::logic_error("Could not find matching feasible stereopermutation for stereocenter");
    }

    mol.assignStereopermutator(permutator.placement(), assignmentIter - std::begin(assignables));
  }
}

void MoleculeBuilder::setBondStereo(
  std::vector<Molecule>& molecules,
  const std::vector<unsigned>& componentMap,
  const std::vector<PrivateGraph::Vertex>& indexInComponentMap,
  const std::string& smiles
) {
  /* Setting the bond stereo from the forward and backward markers is tricky
   * for several reasons.
   *
   * 1. The "/" and "\" markers indicate the "up" or "down" positioning
   * relative to the carbon atom, which may freely occur before or after the
   * marked atom.
   *
   * 2. All we have from the parsing is the sequence in which the
   * stereomarkers were discovered. During parsing, none of these are
   * directly interpreted. We have to make sure they make sense in the first
   * place (no two markers on one side of the double bond may indicate the
   * same relative positioning).
   *
   * That said, I think there are some patterns we can exploit to structure
   * the mess.
   *
   * F/C=C/F is trans, and so is C(\F)=C/F
   *
   * The only side for which there is freedom of reordering for the
   * stereomarkers is left of the bond. I.e. F/C=(\F)C is invalid. So the only
   * side of the double bond for which we have to figure out the relative
   * ordering is left. We can exploit the double bond that must be present in
   * the graph:
   *
   * F/C=C/F -> AB forward, CD forward and BC double bonded
   * C(\F)=C/F -> AB backward, CD forward and AC double bonded
   *
   * It is legal to mark the second substituent at either side too: We can
   * express this molecule also as:
   *
   * [H]\C(\F)=C/F(\[H]) -> AB backward, BC backward, DE forward, DF backward, B = D
   * C(\F)(/[H])=C/F(\[H]) -> AB backward, AC forward, DE forward, DF backward, A = D
   *
   * Unfortunately, all we have is the ordered sequence of marked bonds, so
   * we have to write a state machine that can deal with that. Index
   * repetitions and changes can be our guides to deciding when we have
   * crossed sides of the bond.
   */
  using Temple::Functor::first;
  using Temple::Functor::second;
  auto marker = Temple::Functor::get<2>();

  using Iterator = std::vector<StereoMarkedBondTuple>::const_iterator;
  auto iter = std::cbegin(stereoMarkedBonds);
  const auto end = std::cend(stereoMarkedBonds);

  while(iter != end) {
    SmilesBondStereo state;

    const PrivateGraph::Vertex A = Temple::Functor::first(*iter);
    const PrivateGraph::Vertex B = Temple::Functor::second(*iter);

    // We assume that all vertices are in the same component
    Molecule& mol = molecules.at(componentMap.at(A));

    std::vector<Iterator> leftMarkers {iter};
    std::vector<Iterator> rightMarkers;

    /* Iff the first two marked bonds have an overlapping atom index, then
     * they are on the same side of the bond. The overlapping bond must
     * then be the left atom.
     */
    auto explorer = iter + 1;
    if(explorer == end) {
      throw std::runtime_error("Missing right side of stereo-marked double bond");
    }

    // Check for second marker left of bond
    {
      const PrivateGraph::Vertex X = first(*explorer);

      if(A == X) {
        // Two markers left of bond, C(\F)(/[H]) pattern
        state.left = A;
        leftMarkers.push_back(explorer);
        ++explorer;
      } else if(B == X) {
        // Two markers left of bond, [H]\C(\F) pattern
        state.left = B;
        leftMarkers.push_back(explorer);
        ++explorer;
      }
    }

    // Now we have ensured the explorer is right of the bond
    if(explorer == end) {
      throw std::runtime_error("Missing right side of stereo-marked double bond");
    }

    const auto bondTypeOption = [&](const PrivateGraph::Vertex a, const PrivateGraph::Vertex b) {
      return Temple::Optionals::map(
        mol.graph().bond(
          indexInComponentMap.at(a),
          indexInComponentMap.at(b)
        ),
        [&](const BondIndex& bond) -> BondType {
          return mol.graph().bondType(bond);
        }
      );
    };

    // Establish the right atom
    {
      rightMarkers.push_back(explorer);
      state.right = first(*explorer);

      // Establish the left atom if it is unknown
      if(!state.left) {
        if(bondTypeOption(A, state.right) == BondType::Double) {
          state.left = A;
        } else if(bondTypeOption(B, state.right) == BondType::Double) {
          state.left = B;
        } else {
          throw std::runtime_error("Right side of marked double bond expected, got unrelated bond");
        }
      }
    }

    /* Check for an additional right marker and place explorer at end of
     * relevant marked bonds
     */
    ++explorer;
    if(explorer != end && first(*explorer) == state.right) {
      rightMarkers.push_back(explorer);
      ++explorer;
    }

    assert(state.left);
    /* Now process the collected markers for directionality */
    for(const Iterator& leftMarker : leftMarkers) {
      /* Four cases:
       * - first of the marker is left and forward: C(/F) -> second is up
       * - first of the marker is left and backward: C(\F) -> second is down
       * - second of the marker is left and forward: F/C -> first is down
       * - second of the marker is left and backward: F\C -> first is up
       */
      const bool firstIsLeft = (first(*leftMarker) == state.left.value());
      const bool markerIsForward = (marker(*leftMarker) == SmilesBondType::Forward);
      const bool up = (firstIsLeft == markerIsForward);
      const PrivateGraph::Vertex which = (firstIsLeft ? second(*leftMarker) : first(*leftMarker));

      if(up) {
        if(state.upOfLeft) {
          throw std::runtime_error("Both markers left of double bond indicate 'up' directionality");
        }
        state.upOfLeft = which;
      } else {
        if(state.downOfLeft) {
          throw std::runtime_error("Both markers left of double bond indicate 'down' directionality");
        }
        state.downOfLeft = which;
      }
    }
    for(const Iterator& rightMarker : rightMarkers) {
      assert(first(*rightMarker) == state.right);
      if(marker(*rightMarker) == SmilesBondType::Forward) {
        if(state.upOfRight) {
          throw std::runtime_error("Both markers right of double bond indicate 'up' directionality");
        }
        state.upOfRight = second(*rightMarker);
      } else { // Backward
        if(state.downOfRight) {
          throw std::runtime_error("Both markers right of double bond indicate 'down' directionality");
        }
        state.downOfRight = second(*rightMarker);
      }
    }

    /* Add the information to the molecular graph */
    const auto molBondOption = mol.graph().bond(
      indexInComponentMap.at(state.left.value()),
      indexInComponentMap.at(state.right)
    );
    assert(molBondOption);

    if(auto stereopermutatorOption = mol.stereopermutators().option(molBondOption.value())) {
      if(stereopermutatorOption->numAssignments() == 2) {
        mol.assignStereopermutator(
          molBondOption.value(),
          state.findAssignment(
            *stereopermutatorOption,
            mol,
            indexInComponentMap
          )
        );
      } else {
        std::cerr << "Warning: Smiles '" << smiles << "' contains stereo markers for non-stereogenic double bond\n";
      }
    } else {
      std::cerr << "Warning: Smiles '" << smiles << "' contains stereo markers for non-stereogenic double bond\n";
    }

    // Advance the iterator
    iter = explorer;
  }
}

std::vector<Molecule> MoleculeBuilder::interpret(const std::string& smiles) {
  if(!ringClosures.empty()) {
    throw std::runtime_error("Unmatched ring closure markers remain!");
  }

  std::vector<unsigned> componentMap;
  const unsigned M = graph.connectedComponents(componentMap);

  std::vector<PrivateGraph> precursors;
  precursors.resize(M);

  const unsigned N = graph.V();

  std::vector<PrivateGraph::Vertex> indexInComponentMap(N);
  // Copy vertices
  for(unsigned i = 0; i < N; ++i) {
    auto& precursor = precursors.at(componentMap.at(i));
    PrivateGraph::Vertex newIndex = precursor.addVertex(graph.elementType(i));
    indexInComponentMap.at(i) = newIndex;
  }

  /* Copy edges into the separate components */
  for(const PrivateGraph::Edge& edge : graph.edges()) {
    const PrivateGraph::Vertex source = graph.source(edge);
    const PrivateGraph::Vertex target = graph.target(edge);

    // Both vertices must be in the same component
    auto& precursor = precursors.at(componentMap.at(source));

    precursor.addEdge(
      indexInComponentMap.at(source),
      indexInComponentMap.at(target),
      graph.bondType(edge)
    );
  }

  // Mark eta bonds
  for(auto& precursor : precursors) {
    GraphAlgorithms::updateEtaBonds(precursor);
  }

  const auto aromaticityIncrements = matchAromatics(precursors, componentMap, indexInComponentMap);

  /* Valence fill organic subset in each precursor */
  assert(vertexData.size() == N);
  for(unsigned i = 0; i < N; ++i) {
    const AtomData& data = vertexData.at(i);
    auto& precursor = precursors.at(componentMap.at(i));
    PrivateGraph::Vertex vertexInPrecursor = indexInComponentMap.at(i);

    if(data.hCount) {
      // Fill with specified number of hydrogen atoms
      for(unsigned j = 0; j < data.hCount.value(); ++j) {
        PrivateGraph::Vertex newHydrogenVertex = precursor.addVertex(Utils::ElementType::H);
        precursor.addEdge(vertexInPrecursor, newHydrogenVertex, BondType::Single);
      }
    } else if(!data.atomBracket && isValenceFillElement(precursor.elementType(vertexInPrecursor))) {
      const int currentValence = vertexValence(vertexInPrecursor, precursor);
      const int adjustedValence = aromaticityIncrements.count(i) > 0 ? currentValence + 1 : currentValence;
      const unsigned fillCount = valenceFillElementImplicitHydrogenCount(
        adjustedValence,
        precursor.elementType(vertexInPrecursor)
      );

      for(unsigned j = 0; j < fillCount; ++j) {
        const PrivateGraph::Vertex newHydrogenVertex = precursor.addVertex(Utils::ElementType::H);
        precursor.addEdge(vertexInPrecursor, newHydrogenVertex, BondType::Single);
      }
    }
  }

  /* Convert the graphs to molecules */
  std::vector<Molecule> molecules;
  molecules.reserve(M);
  for(auto&& precursor : precursors) {
    molecules.emplace_back(
      Graph(std::move(precursor))
    );
  }

  /* Stereo routines */
  setShapes(molecules, componentMap, indexInComponentMap);
  setAtomStereo(molecules, componentMap, indexInComponentMap, smiles);
  setBondStereo(molecules, componentMap, indexInComponentMap, smiles);
  addAromaticBondStereo(molecules, componentMap, indexInComponentMap);

  return molecules;
}

std::unordered_set<PrivateGraph::Vertex> MoleculeBuilder::matchAromatics(
  std::vector<PrivateGraph>& precursors,
  const std::vector<unsigned>& componentMap,
  const std::vector<PrivateGraph::Vertex>& indexInComponentMap
) const {
  std::unordered_set<PrivateGraph::Vertex> matchedVertices;

  const unsigned P = precursors.size();
  std::vector<
    std::vector<unsigned>
  > invertedComponentMap(P);
  const unsigned N = vertexData.size();
  for(unsigned i = 0; i < N; ++i) {
    invertedComponentMap.at(componentMap.at(i)).push_back(i);
  }

  // Group piSubgraphEdges by the component they're in
  std::vector<std::vector<PrivateGraph::Edge>> componentPiSubgraphEdges(P);
  for(const auto& edge : piSubgraphEdges) {
    const unsigned componentIndex = componentMap.at(graph.source(edge));
    const auto& precursor = precursors.at(componentIndex);
    componentPiSubgraphEdges.at(componentIndex).push_back(
      precursor.edge(
        indexInComponentMap.at(graph.source(edge)),
        indexInComponentMap.at(graph.target(edge))
      )
    );
  }

  // Construct a pi subgraph for each connected component precursor
  for(unsigned p = 0; p < P; ++p) {
    auto& precursor = precursors.at(p);
    PiSubgraph subgraph;

    /* We're going to iterate over edges, so we'll memoize the viability
     * omissibility check for vertices to avoid repeating work.
     */

    std::unordered_map<PrivateGraph::Vertex, PiSubgraph::EligibleOmissible> viabilityOmissibilityCache;
    const auto memoizedViableOmissible = [&](const PrivateGraph::Vertex i) -> PiSubgraph::EligibleOmissible {
      const auto cached = viabilityOmissibilityCache.find(i);
      if(cached != viabilityOmissibilityCache.end()) {
        return cached->second;
      }

      const PrivateGraph::Vertex j = invertedComponentMap.at(p).at(i);
      const AtomData& atomData = vertexData.at(j);

      const auto viableOmissible = PiSubgraph::eligibleOmissible(i, precursor, atomData);
      viabilityOmissibilityCache.emplace(i, viableOmissible);
      return viableOmissible;
    };

    // Make a pi subgraph with all eligible edges
    for(const auto& edge : componentPiSubgraphEdges.at(p)) {
      const PrivateGraph::Vertex i = precursor.source(edge);
      const PrivateGraph::Vertex j = precursor.target(edge);

      const auto iViableOmissible = memoizedViableOmissible(i);
      const auto jViableOmissible = memoizedViableOmissible(j);

      if(!iViableOmissible.eligible || !jViableOmissible.eligible) {
        continue;
      }

      const auto a = subgraph.findOrAdd(i);
      const auto b = subgraph.findOrAdd(j);

      boost::add_edge(a, b, subgraph.graph);

      if(iViableOmissible.omissible) {
        subgraph.omissible.insert(a);
      }
      if(jViableOmissible.omissible) {
        subgraph.omissible.insert(b);
      }
    }

    if(boost::num_vertices(subgraph.graph) > 0) {
      if(const auto subgraphMatches = subgraph.match()) {
        for(PrivateGraph::Vertex i : subgraphMatches.value()) {
          matchedVertices.insert(invertedComponentMap.at(p).at(i));
        }
      } else {
        /* Matching must succeed for each component's pi subgraph, otherwise
         * the smiles is invalid.
         */
        throw std::runtime_error("Aromatic smiles subgraph could not be matched!");
      }
    }
  }

  return matchedVertices;
}

void MoleculeBuilder::addAromaticBondStereo(
  std::vector<Molecule>& molecules,
  const std::vector<unsigned>& componentMap,
  const std::vector<PrivateGraph::Vertex>& indexInComponentMap
) {
  std::vector<std::unordered_set<AtomIndex>> aromaticAtoms(molecules.size());
  for(AtomIndex i = 0; i < vertexData.size(); ++i) {
    const AtomData& data = vertexData.at(i);
    if(data.partialElement.aromatic) {
      aromaticAtoms.at(componentMap.at(i)).insert(indexInComponentMap.at(i));
    }
  }

  for(unsigned i = 0; i < molecules.size(); ++i) {
    Molecule& mol = molecules.at(i);
    const auto& aromatics = aromaticAtoms.at(i);

    for(const auto& cycleEdges : mol.graph().cycles()) {
      const bool allAromatic = Temple::all_of(
        cycleEdges,
        [&](const BondIndex edge) -> bool {
          return aromatics.count(edge.first) > 0;
        }
      );

      if(!allAromatic) {
        continue;
      }

      for(const auto& bond : cycleEdges) {
        auto permutator = mol.stereopermutators().option(bond);
        if(!permutator) {
          mol.addPermutator(bond);
        }
      }
    }
  }
}

} // namespace IO
} // namespace Molassembler
} // namespace Scine
