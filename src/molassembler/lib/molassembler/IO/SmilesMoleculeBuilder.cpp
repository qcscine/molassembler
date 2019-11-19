#include "molassembler/IO/SmilesMoleculeBuilder.h"

#include "molassembler/StereopermutatorList.h"
#include "molassembler/Stereopermutators/AbstractPermutations.h"
#include "molassembler/Stereopermutators/FeasiblePermutations.h"
#include "molassembler/Stereopermutators/ShapeVertexMaps.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Modeling/BondDistance.h"
#include "molassembler/Modeling/LocalGeometryModel.h"
#include "molassembler/Molecule.h"
#include "shapes/Data.h"

#include "temple/Optionals.h"
#include "temple/Functional.h"

#include <iostream>

namespace Scine {
namespace molassembler {
namespace IO {

bool MoleculeBuilder::isValenceFillElement(Utils::ElementType e) {
  const unsigned Z = Utils::ElementInfo::Z(e);
  if(5 <= Z && Z <= 9) {
    // B, C, N, O, F
    return true;
  }

  if(15 <= Z && Z <= 17) {
    // P, S, Cl
    return true;
  }

  if(Z == 35 || Z == 53) {
    // Br, I
    return true;
  }

  return false;
}

unsigned MoleculeBuilder::valenceFillElementImplicitHydrogenCount(
  const int valence,
  Utils::ElementType e
) {
  assert(valence >= 0);
  assert(isValenceFillElement(e));

  /* Quoting from the spec:
   *
   * The implicit hydrogen count is determined by summing the bond orders of
   * the bonds connected to the atom. If that sum is equal to a known valence
   * for the element or is greater than any known valence then the implicit
   * hydrogen count is 0. Otherwise the implicit hydrogen count is the
   * difference between that sum and the next highest known valence.
   */

  switch(Utils::ElementInfo::Z(e)) {
    case 5: return std::max(0, 3 - valence); // B
    case 6: return std::max(0, 4 - valence); // C
    case 7: { // N
      return std::min(
        std::max(0, 3 - valence),
        std::max(0, 5 - valence)
      );
    }
    case 8: return std::max(0, 2 - valence); // O
    case 15: { // P
      return std::min(
        std::max(0, 3 - valence),
        std::max(0, 5 - valence)
      );
    }
    case 16: { // S
      return std::min({
        std::max(0, 2 - valence),
        std::max(0, 4 - valence),
        std::max(0, 6 - valence)
      });
    }
    default: return std::max(0, 1 - valence); // F, Cl, Br, I are the remaining cases
  }
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

std::vector<unsigned> MoleculeBuilder::shapeMap(const ChiralData& chiralData) {
  if(chiralData.shape == Shapes::Shape::Tetrahedron) {
    switch(chiralData.chiralIndex) {
      case 1: return {{0, 1, 2, 3}}; // @, TH1
      case 2: return {{0, 1, 3, 2}}; // @@, TH2
    }
  } else if(chiralData.shape == Shapes::Shape::Square) {
    switch(chiralData.chiralIndex) {
      case 1: return {{0, 1, 2, 3}}; // SP1 = U
      case 2: return {{0, 2, 3, 1}}; // SP2 = 4
      case 3: return {{3, 2, 0, 1}}; // SP3 = Z
    }
  } else if(chiralData.shape == Shapes::Shape::TrigonalBipyramid) {
    switch(chiralData.chiralIndex) {
      case  1: return {{1, 2, 3, 0, 4}}; // TB1 = a, e, @
      case  2: return {{1, 3, 2, 0, 4}}; // TB2 = a, e, @@
      case  3: return {{1, 2, 4, 0, 3}}; // TB3 = a, d, @
      case  4: return {{1, 4, 2, 0, 3}}; // TB4 = a, d, @@
      case  5: return {{1, 3, 4, 0, 2}}; // TB5 = a, c, @
      case  6: return {{1, 4, 3, 0, 2}}; // TB6 = a, c, @@
      case  7: return {{2, 3, 4, 0, 1}}; // TB7 = a, b, @
      case  8: return {{2, 4, 3, 0, 1}}; // TB8 = a, b, @@
      case  9: return {{0, 2, 3, 1, 4}}; // TB9 = b, e, @
      case 10: return {{0, 2, 4, 1, 3}}; // TB10 = b, d, @
      case 11: return {{0, 3, 2, 1, 4}}; // TB11 = b, e, @@
      case 12: return {{0, 4, 2, 1, 3}}; // TB12 = b, d, @@
      case 13: return {{0, 3, 4, 1, 2}}; // TB13 = b, c, @
      case 14: return {{0, 4, 3, 1, 2}}; // TB14 = b, c, @@
      case 15: return {{0, 1, 3, 2, 4}}; // TB15 = c, e, @
      case 16: return {{0, 1, 4, 2, 3}}; // TB16 = c, d, @
      case 17: return {{0, 1, 2, 3, 4}}; // TB17 = d, e, @
      case 18: return {{0, 2, 1, 3, 4}}; // TB18 = d, e, @@
      case 19: return {{0, 4, 1, 2, 3}}; // TB19 = c, d, @@
      case 20: return {{0, 3, 1, 2, 4}}; // TB20 = c, e, @@
    }
  } else if(chiralData.shape == Shapes::Shape::Octahedron) {
    switch(chiralData.chiralIndex) {
      /* Look along an axis, what remains is a square with a winding. So
       * square shapes are reused with definitions of square shapes and windings
       */
      case  1: return {{1, 2, 3, 4, 0, 5}}; // OH1 = a, f, U, @
      case  2: return {{4, 3, 2, 1, 0, 5}}; // OH2 = a, f, U, @@
      case  3: return {{1, 2, 3, 5, 0, 4}}; // OH3 = a, e, U, @
      case 16: return {{5, 3, 2, 1, 0, 4}}; // OH16 = a, e, U, @@
      case  6: return {{1, 2, 4, 5, 0, 3}}; // OH6 = a, d, U, @
      case 18: return {{5, 4, 2, 1, 0, 3}}; // OH18 = a, d, U, @@
      case 19: return {{1, 3, 4, 5, 0, 2}}; // OH19 = a, c, U, @
      case 24: return {{5, 4, 3, 1, 0, 2}}; // OH24 = a, c, U, @@
      case 25: return {{2, 3, 4, 5, 0, 1}}; // OH25 = a, b, U, @
      case 30: return {{5, 4, 3, 2, 0, 1}}; // OH30 = a, b, U, @@

      /* Note: For the Z shape, the connection between the first two atoms
       * determines the winding.
       */
      case  4: return {{1, 2, 4, 3, 0, 5}}; // OH4 = a, f, Z, @
      case 14: return {{3, 4, 2, 1, 0, 5}}; // OH14 = a, f, Z, @@
      case  5: return {{1, 2, 5, 3, 0, 4}}; // OH5 = a, e, Z, @
      case 15: return {{3, 5, 2, 1, 0, 4}}; // OH15 = a, e, Z, @@
      case  7: return {{1, 2, 5, 4, 0, 3}}; // OH7 = a, d, Z, @
      case 17: return {{4, 5, 2, 1, 0, 3}}; // OH17 = a, d, Z, @@
      case 20: return {{1, 3, 5, 4, 0, 2}}; // OH20 = a, c, Z, @
      case 23: return {{4, 5, 3, 1, 0, 2}}; // OH23 = a, c, Z, @@
      case 26: return {{2, 3, 5, 4, 0, 1}}; // OH26 = a, b, Z, @
      case 29: return {{4, 5, 3, 2, 0, 1}}; // OH29 = a, b, Z, @@

      /* For the 4 shape, the connection between the second and third atom
       * determines the winding.
       */
      case 10: return {{4, 2, 3, 1, 0, 5}}; // OH10 = a, f, 4, @
      case  8: return {{1, 3, 2, 4, 0, 5}}; // OH8 = a, f, 4, @@
      case 11: return {{5, 2, 3, 1, 0, 4}}; // OH11 = a, e, 4, @
      case  9: return {{1, 3, 2, 5, 0, 4}}; // OH9 = a, e, 4, @@
      case 13: return {{5, 2, 4, 1, 0, 3}}; // OH13 = a, d, 4, @
      case 12: return {{1, 4, 2, 5, 0, 3}}; // OH12 = a, d, 4, @@
      case 22: return {{5, 3, 4, 1, 0, 2}}; // OH22 = a, c, 4, @
      case 21: return {{1, 4, 3, 5, 0, 2}}; // OH21 = a, c, 4, @@
      case 28: return {{5, 3, 4, 2, 0, 1}}; // OH28 = a, b, 4, @
      case 27: return {{2, 4, 3, 5, 0, 1}}; // OH27 = a, b, 4, @@
    }
  }

  throw std::logic_error("Invalid combination of shape and chiral index!");
}

void MoleculeBuilder::addAtom(const AtomData& atom) {
  InnerGraph::Vertex newVertex = graph.addVertex(atom.getElement());

  if(atom.partialElement.Z == 1 && atom.hCount && atom.hCount.value() != 0) {
    throw std::runtime_error("Hydrogen atoms cannot have hydrogen counts!");
  }

  vertexData.push_back(atom);

  if(lastBondData.which() == 0) {
    auto data = boost::get<SimpleLastBondData>(lastBondData);
    if(data == SimpleLastBondData::Unspecified) {
      assert(!vertexStack.empty());
      graph.addEdge(
        vertexStack.top(),
        newVertex,
        BondType::Single
      );
    }
  } else {
    assert(!vertexStack.empty());
    auto data = boost::get<BondData>(lastBondData);
    graph.addEdge(
      vertexStack.top(),
      newVertex,
      data.type.value_or(BondType::Single)
    );

    // Store stereo-marked bonds for later
    if(data.ezStereo) {
      stereoMarkedBonds.emplace_back(
        vertexStack.top(),
        newVertex,
        data.ezStereo.value()
      );
    }
  }

  if(vertexStack.empty()) {
    vertexStack.push(newVertex);
  } else {
    vertexStack.top() = newVertex;
  }

  lastBondData = SimpleLastBondData::Unspecified;
}

void MoleculeBuilder::addRingClosure(const BondData& bond) {
  assert(bond.ringNumber);
  const unsigned key = bond.ringNumber.value();
  auto findIter = ringClosures.find(key);
  if(findIter == std::end(ringClosures)) {
    // Add the entry to the map for later
    ringClosures.emplace_hint(
      findIter,
      std::piecewise_construct,
      std::make_tuple(key),
      std::make_tuple(vertexStack.top(), bond.type)
    );
  } else {
    // Add the edge to the graph now and remove the map entry
    InnerGraph::Vertex a = findIter->second.first;
    InnerGraph::Vertex b = vertexStack.top();

    if(a == b) {
      throw std::runtime_error("Same-atom ring-closing bond!");
    }

    if(graph.edgeOption(a, b) != boost::none) {
      throw std::runtime_error("Ring closing bond already exists!");
    }

    // Ensure the specified bond types match (this fn can throw)
    const BondType type = mutualBondType(
      findIter->second.second,
      bond.type
    );

    graph.addEdge(a, b, type);

    // Remove the entry from the map
    ringClosures.erase(findIter);
  }
}

void MoleculeBuilder::setShapes(
  std::vector<Molecule>& molecules,
  const std::vector<unsigned>& componentMap,
  const std::vector<InnerGraph::Vertex>& indexInComponentMap
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

      auto modelArgs = LocalGeometry::reduceToSiteInformation(
        mol.graph(),
        atomIndex,
        stereoOption->getRanking()
      );

      auto shapeOption = LocalGeometry::vsepr(
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
  const std::vector<InnerGraph::Vertex>& indexInComponentMap
) {
  const unsigned N = vertexData.size();
  for(unsigned i = 0; i < N; ++i) {
    const AtomData& atomData = vertexData.at(i);
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
      std::cerr << "Warning: Smiles contains a stereo marker for a non-stereogenic " << Shapes::name(chiralData.shape) << " shape center\n";
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
    std::vector<unsigned> sortedSites = temple::sort(
      temple::iota<unsigned>(S),
      [&](const unsigned a, const unsigned b) -> bool {
        return ranking.sites.at(a) < ranking.sites.at(b);
      }
    );

    /* Atom bracket hcount special case:
     *
     * If one of the neighbors is a hydrogen atom and is represented as an
     * hcount instead of explicitly, then it is considered to be the first
     * atom in the clockwise or anticlockwise counting.
     */
    for(unsigned j = 0; j < ranking.sites.size(); ++j) {
      if(
        ranking.sites.at(j).size() == 1
        && mol.graph().elementType(ranking.sites.at(j).front()) == Utils::ElementType::H
        && vertexData.at(i).hCount == 1u
      ) {
        auto vertexMapIter = std::find(std::begin(sortedSites), std::end(sortedSites), j);
        assert(vertexMapIter != std::end(sortedSites));
        std::iter_swap(std::begin(sortedSites), vertexMapIter);
        break;
      }
    }

    // TODO missing weirdness re: ring closing bonds

    /* Transfer the sorted sites onto shape vertices */
    auto vertexMap = shapeMap(chiralData);
    std::vector<unsigned> siteToShapeVertexMap(S);
    for(unsigned j = 0; j < S; ++j) {
      siteToShapeVertexMap.at(vertexMap.at(j)) = sortedSites.at(j);
      // siteToShapeVertexMap.at(j) = sortedSites.at(vertexMap.at(j));
    }

    /* Create a stereopermutation and look for it in the list of feasibles */
    auto soughtStereopermutation = stereopermutationFromSiteToShapeVertexMap(
      siteToShapeVertexMap,
      ranking.links,
      permutator.getAbstract().canonicalSites
    );

    auto soughtRotations = soughtStereopermutation.generateAllRotations(chiralData.shape);
    const auto& assignables = permutator.getFeasible().indices;
    auto assignmentIter = temple::find_if(
      assignables,
      [&](const unsigned stereopermutationIndex) -> bool {
        const auto& stereopermutation = permutator.getAbstract().permutations.stereopermutations.at(stereopermutationIndex);
        return soughtRotations.count(stereopermutation) > 0;
      }
    );

    if(assignmentIter == std::end(assignables)) {
      throw std::logic_error("Could not find matching feasible stereopermutation for stereocenter");
    }

    mol.assignStereopermutator(permutator.centralIndex(), assignmentIter - std::begin(assignables));
  }
}

void MoleculeBuilder::setBondStereo(
  std::vector<Molecule>& molecules,
  const std::vector<unsigned>& componentMap,
  const std::vector<InnerGraph::Vertex>& indexInComponentMap
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
  struct BondStereo {
    boost::optional<InnerGraph::Vertex> left;
    InnerGraph::Vertex right;
    boost::optional<InnerGraph::Vertex> upOfLeft;
    boost::optional<InnerGraph::Vertex> downOfLeft;
    boost::optional<InnerGraph::Vertex> upOfRight;
    boost::optional<InnerGraph::Vertex> downOfRight;

    unsigned findAssignment(
      BondStereopermutator stereopermutator,
      const Molecule& mol,
      const std::vector<InnerGraph::Vertex>& indexInComponentMap
    ) const {
      auto first = mol.stereopermutators().option(stereopermutator.edge().first).value();
      auto second = mol.stereopermutators().option(stereopermutator.edge().second).value();

      if(first.centralIndex() == indexInComponentMap.at(right)) {
        std::swap(first, second);
      }

      auto getSiteIndexLeft = [&](const InnerGraph::Vertex i) -> unsigned {
        return first.getRanking().getSiteIndexOf(indexInComponentMap.at(i));
      };
      auto getSiteIndexRight = [&](const InnerGraph::Vertex i) -> unsigned {
        return second.getRanking().getSiteIndexOf(indexInComponentMap.at(i));
      };

      assert(stereopermutator.numAssignments() == 2);
      for(unsigned i = 0; i < 2; ++i) {
        stereopermutator.assign(i);

        auto upOfLeftSiteIndex = temple::optionals::map(upOfLeft, getSiteIndexLeft);
        auto downOfLeftSiteIndex = temple::optionals::map(downOfLeft, getSiteIndexLeft);
        auto upOfRightSiteIndex = temple::optionals::map(upOfRight, getSiteIndexRight);
        auto downOfRightSiteIndex = temple::optionals::map(downOfRight, getSiteIndexRight);

        if(upOfLeftSiteIndex) {
          if(upOfRightSiteIndex) {
            if(std::fabs(stereopermutator.dihedral(first, *upOfLeftSiteIndex, second, *upOfRightSiteIndex)) > 1e-3) {
              continue;
            }
          }
          if(downOfRightSiteIndex) {
            if(std::fabs(stereopermutator.dihedral(first, *upOfLeftSiteIndex, second, *downOfRightSiteIndex) - M_PI) > 1e-3) {
              continue;
            }
          }
        }
        if(downOfLeftSiteIndex) {
          if(upOfRightSiteIndex) {
            if(std::fabs(stereopermutator.dihedral(first, *downOfLeftSiteIndex, second, *upOfRightSiteIndex) - M_PI) > 1e-3) {
              continue;
            }
          }
          if(downOfRightSiteIndex) {
            if(std::fabs(stereopermutator.dihedral(first, *downOfLeftSiteIndex, second, *downOfRightSiteIndex)) > 1e-3) {
              continue;
            }
          }
        }

        return i;
      }

      throw std::logic_error("Failed to find matching stereopermutation for BondStereo state.");
    }
  };

  auto first = [](const auto& tup) { return std::get<0>(tup); };
  auto second = [](const auto& tup) { return std::get<1>(tup); };
  auto marker = [](const auto& tup) { return std::get<2>(tup); };

  using Iterator = std::vector<StereoMarkedBondTuple>::const_iterator;
  Iterator start = std::cbegin(stereoMarkedBonds);
  const Iterator end = std::cend(stereoMarkedBonds);

  while(start != end) {
    BondStereo state;

    InnerGraph::Vertex A = first(*start);
    InnerGraph::Vertex B = second(*start);

    // We assume that all vertices are in the same component
    Molecule& mol = molecules.at(componentMap.at(A));

    std::vector<Iterator> leftMarkers {start};
    std::vector<Iterator> rightMarkers;

    /* Iff the first two marked bonds have an overlapping atom index, then
     * they are on the same side of the bond. The overlapping bond must
     * then be the left atom.
     */
    Iterator explorer = start + 1;
    if(explorer == end) {
      throw std::runtime_error("Missing right side of stereo-marked double bond");
    }

    // Check for second marker left of bond
    {
      InnerGraph::Vertex X = first(*explorer);

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

    auto bondTypeOption = [&](const InnerGraph::Vertex a, const InnerGraph::Vertex b) {
      return temple::optionals::map(
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
      bool firstIsLeft = (first(*leftMarker) == state.left.value());
      bool markerIsForward = (marker(*leftMarker) == BondData::StereoMarker::Forward);
      bool up = (firstIsLeft == markerIsForward);
      InnerGraph::Vertex which = (firstIsLeft ? second(*leftMarker) : first(*leftMarker));

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
      if(marker(*rightMarker) == BondData::StereoMarker::Forward) {
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
    auto molBondOption = mol.graph().bond(
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
        std::cerr << "Warning: Smiles contains stereo markers for non-stereogenic double bond\n";
      }
    } else {
      std::cerr << "Warning: Smiles contains stereo markers for non-stereogenic double bond\n";
    }

    // Advance the start iterator
    start = explorer;
  }
}

std::vector<Molecule> MoleculeBuilder::interpret() {
  if(!ringClosures.empty()) {
    throw std::runtime_error("Unmatched ring closure markers remain!");
  }

  std::vector<unsigned> componentMap;
  const unsigned M = graph.connectedComponents(componentMap);

  std::vector<InnerGraph> precursors;
  precursors.resize(M);

  const unsigned N = graph.N();

  std::vector<InnerGraph::Vertex> indexInComponentMap(N);
  // Copy vertices
  for(unsigned i = 0; i < N; ++i) {
    auto& precursor = precursors.at(componentMap.at(i));
    InnerGraph::Vertex newIndex = precursor.addVertex(graph.elementType(i));
    indexInComponentMap.at(i) = newIndex;
  }

  /* Copy edges into the separate components */
  for(const InnerGraph::Edge& edge : boost::make_iterator_range(graph.edges())) {
    const InnerGraph::Vertex source = graph.source(edge);
    const InnerGraph::Vertex target = graph.target(edge);

    // Both vertices must be in the same component
    auto& precursor = precursors.at(componentMap.at(source));

    precursor.addEdge(
      indexInComponentMap.at(source),
      indexInComponentMap.at(target),
      graph.bondType(edge)
    );
  }

  /* Valence fill organic subset in each precursor */
  assert(vertexData.size() == N);
  for(unsigned i = 0; i < N; ++i) {
    const AtomData& data = vertexData.at(i);
    auto& precursor = precursors.at(componentMap.at(i));
    InnerGraph::Vertex vertexInPrecursor = indexInComponentMap.at(i);

    if(data.hCount) {
      // Fill with specified number of hydrogen atoms
      for(unsigned j = 0; j < data.hCount.value(); ++j) {
        InnerGraph::Vertex newHydrogenVertex = precursor.addVertex(Utils::ElementType::H);
        precursor.addEdge(vertexInPrecursor, newHydrogenVertex, BondType::Single);
      }
    } else if(!data.atomBracket && isValenceFillElement(precursor.elementType(vertexInPrecursor))) {
      // Figure out current valence.
      int currentValence = 0;
      for(
        const InnerGraph::Edge edge :
        boost::make_iterator_range(precursor.edges(vertexInPrecursor))
      ) {
        currentValence += Bond::bondOrderMap.at(
          static_cast<unsigned>(
            precursor.bondType(edge)
          )
        );
      }

      const unsigned fillCount = valenceFillElementImplicitHydrogenCount(
        currentValence,
        precursor.elementType(vertexInPrecursor)
      );

      for(unsigned j = 0; j < fillCount; ++j) {
        InnerGraph::Vertex newHydrogenVertex = precursor.addVertex(Utils::ElementType::H);
        precursor.addEdge(vertexInPrecursor, newHydrogenVertex, BondType::Single);
      }
    }
  }

  /* Convert the graphs to molecules */
  std::vector<Molecule> molecules;
  molecules.reserve(M);
  for(auto&& precursor : precursors) {
    molecules.emplace_back(
      OuterGraph(std::move(precursor))
    );
  }

  /* Stereo routines */
  setShapes(molecules, componentMap, indexInComponentMap);
  setAtomStereo(molecules, componentMap, indexInComponentMap);
  setBondStereo(molecules, componentMap, indexInComponentMap);

  return molecules;
}

} // namespace IO
} // namespace molassembler
} // namespace Scine
