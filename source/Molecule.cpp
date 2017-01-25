#include "Molecule.h"
#include "CN4Stereocenter.h"
#include "BondDistance.h"
#include "symmetry_information/Symmetries.h"
#include "CommonTrig.h"

#include "GraphDistanceMatrix.h"
#include "AdjacencyListAlgorithms.h"
#include "AdjacencyMatrix.h"

#include "DistanceGeometry/DistanceBoundsMatrix.h"

// Delib
#include "ElementInfo.h"

/* TODO
 * - implement remaining functions, importantly the operator ==
 * - consider how to interface GraphFeatures with Molecule, there is a need for
 *   information flow between both:: Molecule operator == needs to be able
 *   to compare GraphFeatures and GraphFeatures need access to element types
 *   and edges
 * - Molecule operator == can be implemented in myriad ways, suggest some 
 *   heuristics before going for testing graph isomorphism:
 *
 *   1. #vertices, #edges O(1)
 *   2. comparison of GraphFeatures 
 *      - same amount
 *      - exactly the same
 *   3. isomorphism
 *
 *   maybe also consider the Wiener index, global clustering and algebraic
 *   connectivity
 *
 */

namespace MoleculeManip {

/* Constructors */
Molecule::Molecule(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  // store ElementTypes
  _elements.push_back(a);
  _elements.push_back(b);
  // update _adjacencies
  _adjacencies.addSlot();
  _adjacencies.addSlot();
  _adjacencies.addAdjacency(0, 1);
  // update _edges
  _edges.add(
    0,
    1,
    bondType
  );
}

Molecule::Molecule(
  const Delib::ElementTypeCollection& elements,
  const AdjacencyList& adjacencies,
  const Edges& edges
) : 
  _elements(elements),
  _adjacencies(adjacencies),
  _edges(edges) { 

  _detectStereocenters();
}

Molecule::Molecule(
  const Delib::ElementTypeCollection& elements,
  const Delib::PositionCollection& positions,
  const AdjacencyList& adjacencies,
  const Edges& edges
) : 
  _elements(elements),
  _positions(positions),
  _adjacencies(adjacencies),
  _edges(edges) {
  
  _detectStereocenters();
}

/* Private member functions */
void Molecule::_detectStereocenters() {
  for(unsigned i = 0; i < _adjacencies.size(); i++) {
    if(
      _adjacencies[i].size() == 4 
      && hydrogenCount(i) < 2 // reduces amount of superfluous constructions
    ) {
      std::shared_ptr<
        Stereocenters::CN4Stereocenter
      > newStereocenter = std::make_shared<
        Stereocenters::CN4Stereocenter
      >(
        this,
        i
      );

      if(newStereocenter -> assignments() > 1) {
        _stereocenters.add(
          std::move(newStereocenter)
        );
      }
    }
  }
}

bool Molecule::_validAtomIndex(const AtomIndexType& a) const {
  return (
    a < _adjacencies.size()
    && _elements[a] != Delib::ElementType::none
  );
}

bool Molecule::_validAtomIndices(
  const AtomIndexType& a,
  const AtomIndexType& b
) const {
  return (
    _validAtomIndex(a)
    && _validAtomIndex(b)
    && a < b
  );
}

AtomIndexType Molecule::addAtom(
  const Delib::ElementType& elementType,
  const AtomIndexType& bondedToIndex,
  const BondType& bondType
) {
  auto addedIndex = _adjacencies.addSlot();

  assert(bondedToIndex < addedIndex);
  _adjacencies.addAdjacency(
    bondedToIndex,
    addedIndex
  );

  _edges.add(
    bondedToIndex,
    addedIndex,
    bondType
  );

  _elements.push_back(elementType);

  return addedIndex;
}

void Molecule::addBond(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const BondType& bondType
) {
  assert(_validAtomIndices(a, b));
  _adjacencies.addAdjacency(
    a,
    b
  );

  _edges.add(
    a,
    b,
    bondType
  );
}

void Molecule::removeAtom(const AtomIndexType& a) {
  assert(_validAtomIndex(a));

  // must update edges and adjacencies 
  auto bonded_to = _adjacencies.getAdjacencies(a);

  // erase all edges to and from this atom
  for(const auto& bondedAtomIndex : bonded_to) {
    _edges.remove(
      std::min(a, bondedAtomIndex),
      std::max(a, bondedAtomIndex)
    );
  }

  // remove all other mentions in _adjacencies
  for(const auto& bondedAtomIndex : bonded_to) {
    _adjacencies.removeAdjacency(a, bondedAtomIndex);
  }

  // set _atom_exists to false for this index
  _elements[a] = Delib::ElementType::none;
}

void Molecule::removeBond(
  const AtomIndexType& a,
  const AtomIndexType& b
) {
  assert(_validAtomIndices(a, b));
  _adjacencies.removeAdjacency(a, b);

  _edges.remove(
    a,
    b
  );
}

Delib::ElementType Molecule::getElementType(
  const AtomIndexType& a
) const {
  assert(_validAtomIndex(a));
  return _elements[a];
}

AtomIndexType Molecule::getNumAtoms() const {
  return _elements.size();
}

EdgeIndexType Molecule::getNumBonds() const {
  return _edges.size();
}

const Edges& Molecule::getEdges() const {
  return _edges;
}

BondType Molecule::getBondType(
  const AtomIndexType& a,
  const AtomIndexType& b
) const {
  auto edgeOption = _edges.get(a, b);
  assert(edgeOption);
  return edgeOption.value();
}

unsigned Molecule::hydrogenCount(const AtomIndexType& a) const {
  assert(_validAtomIndex(a));
  unsigned count = 0;
  for(const auto& adjacentIndex : _adjacencies[a]) {
    if(getElementType(adjacentIndex) == Delib::ElementType::H) {
      count += 1;
    }
  }
  return count;
}

std::vector<AtomIndexType> Molecule::getBondedAtomIndices(
  const AtomIndexType& a
) const {
  assert(_validAtomIndex(a));
  return _adjacencies.getAdjacencies(a);
}

/* TODO
 * - does not treat correctly:
 *   - cycles
 *   - stereocenters (Z over E, R over S (?))
 *   - double and triple bond ghost atom splitting
 * - unsure about sub-lists. is this approach even remotely correct?
 * - FUCK CIP rules -> maybe just use the unsigned values of assignments in
 *   GraphFeatures and rank branches with that.
 * - test
 */
std::pair<
  std::vector<AtomIndexType>, // the sorted list of substituent priorities
  std::set< // a set of pairs of AtomIndexTypes that are EQUAL
    std::pair<
      AtomIndexType,
      AtomIndexType
    >
  >
> Molecule::rankPriority(
  const AtomIndexType& a,
  const std::vector<AtomIndexType>& excludeAdjacent 
) const {
  auto toRank = _adjacencies.getAdjacencies(a);
  std::set<
    std::pair<
      AtomIndexType,
      AtomIndexType
    >
  > equalPairs;

  // remove excludes from toRank
  toRank.erase(
    std::remove_if(
      toRank.begin(),
      toRank.end(),
      [&excludeAdjacent](const AtomIndexType& atomIndex) {
        return std::find(
          excludeAdjacent.begin(),
          excludeAdjacent.end(),
          atomIndex
        ) != excludeAdjacent.end();
      }
    ),
    toRank.end()
  );

  auto getZ = [this](const AtomIndexType& a) -> int {
    return static_cast<int>(getElementType(a));
  };

  auto BFSIterate = [this, getZ](
    std::set<AtomIndexType>& visitedSet,
    std::vector<AtomIndexType>& seeds,
    std::multiset<int, std::greater<int> >& Zs
  ) {
    std::vector<AtomIndexType> newSeeds;
    for(const auto& index : seeds) {
      auto adjacent = getBondedAtomIndices(index);
      for(const auto& potentialSeed : adjacent) {
        // if not in visitedSet
        if(visitedSet.count(potentialSeed) == 0) {
          // add it to the set and new seeds
          visitedSet.insert(potentialSeed);
          newSeeds.push_back(potentialSeed);

          // add it's Z to Zs
          Zs.insert(getZ(potentialSeed));
        } // else skip
      }
    }

    // overwrite seeds
    seeds = newSeeds;
  };
  
  // sort toRank according to CIP-like rules
  std::sort(
    toRank.begin(),
    toRank.end(),
    [&a, this, &getZ, &BFSIterate, &equalPairs](
      const AtomIndexType& lhs,
      const AtomIndexType& rhs
    ) {
      std::set<AtomIndexType> lhsVisited = {a}, rhsVisited = {a};
      std::vector<AtomIndexType> lhsSeeds = {lhs}, rhsSeeds = {rhs};
      std::multiset<
        int,
        std::greater<int> // in CIP, list of Z is ordered DESC
      > lhsZ = { getZ(lhs) }, rhsZ = { getZ(rhs) };
      while(
          lhsSeeds.size() > 0
          || rhsSeeds.size() > 0
      ) {
        // compare lists
        if(lhsZ < rhsZ) return true;
        else if(lhsZ > rhsZ) return false;

        // iterate along the bonds
        BFSIterate(lhsVisited, lhsSeeds, lhsZ);
        BFSIterate(rhsVisited, rhsSeeds, rhsZ);
      }

      // all equal -> add to equalPairs
      equalPairs.emplace(
        std::min(lhs, rhs),
        std::max(lhs, rhs)
      );
      return false;
    }
  );

  return {
    toRank,
    equalPairs
  };
}

void Molecule::_dumpGraphviz(std::ostream& os) const {

  std::map<
    std::string,
    std::string
  > elementBGColorMap {
    {"H", "white"},
    {"C", "gray"},
    {"N", "blue"},
    {"O", "red"}
  };

  std::map<
    std::string,
    std::string
  > elementFGColorMap {
    {"H", "black"},
    {"C", "white"},
    {"N", "white"},
    {"O", "white"}
  };

  auto getSymbolString = [&](const AtomIndexType& index) -> std::string {
    return Delib::ElementInfo::symbol(_elements.at(index));
  };

  auto nodeLabel = [&](const AtomIndexType& index) {
    auto symbol = getSymbolString(index);
    if(elementBGColorMap.count(symbol) > 0) {
      return std::string("\"")
        + std::to_string(index)
        + "\"";
    } else {
      return std::string("\"")
        + std::to_string(index)
        + symbol
        + "\"";
    }
  };

  auto nodeProperties = [&](const std::string& symbolString) -> std::string {
    if(symbolString == "H") {
      return std::string(" [fillcolor=")
        + elementBGColorMap.at("H")
        + ", fontcolor="
        + elementFGColorMap.at("H")
        + ", fontsize=10, width=.3, fixedsize=true]";
    } else if(elementBGColorMap.count(symbolString) == 1) {
      return std::string(" [fillcolor=")
        + elementBGColorMap.at(symbolString)
        + ", fontcolor="
        + elementFGColorMap.at(symbolString)
        + "]";
    } else return " [fillcolor=white, fontcolor=black]";
  };

  auto edgeProperties = [&](
    const std::string& symbolFrom,
    const std::string& symbolTo
  ) {
    if(
      symbolFrom == "H" 
      || symbolTo == "H"
    ) {
      return " [len=0.5]";
    } else {
      return "";
    }
  };

  os << "graph G {\n  graph [fontname = \"Arial\", layout = neato];\n"
    << "  node [fontname = \"Arial\", shape = circle, style = filled];\n"
    << "  edge [fontname = \"Arial\"];\n";

  std::map<
    std::string,
    std::vector<
      std::string
    >
  > elementNameLabelListMap;

  // group elements together
  for(unsigned i = 0; i < _elements.size(); i++) {
    if(elementNameLabelListMap.count(getSymbolString(i)) == 0) {
      elementNameLabelListMap[
        getSymbolString(i)
      ] = {
        nodeLabel(i)
      };
    } else {
      elementNameLabelListMap.at(
        getSymbolString(i)
      ).push_back(
        nodeLabel(i)
      );
    }
  }

  for(const auto& symbolLabelListPair : elementNameLabelListMap) {
    os << "  node" << nodeProperties(symbolLabelListPair.first) << "\n ";
    for(const auto& label : symbolLabelListPair.second) {
      os << " " << label;
    }
    os << ";\n";
  }

  for(const auto& edge: _edges) {
    os << "\n  " << nodeLabel(edge.first.first) << " -- " 
      << nodeLabel(edge.first.second)
      << edgeProperties(
        getSymbolString(edge.first.first),
        getSymbolString(edge.first.second)
      ) 
      << ";";
  }

  os << "\n}\n";
}

std::ostream& operator << (
  std::ostream& os,
  const Molecule& mol
) {
  os << "Begin graphviz –––––––––––––––\n\n";
  mol._dumpGraphviz(os);
  os << "\n––––––––––––––––– End graphviz\n\n";

  if(mol._stereocenters.size() > 0) {
    os << "Stereocenter information:\n";
    for(const auto& stereocenterPtr: mol._stereocenters) {
      os << stereocenterPtr << std::endl;
    }
  }
  
  return os;
}

DistanceGeometry::DistanceBoundsMatrix Molecule::getDistanceBoundsMatrix() const {

  const double oneTwoDelta = 0.01;

  DistanceGeometry::DistanceBoundsMatrix distanceBounds(getNumAtoms());

  /* enter all 1-2 distance bounds into the matrix for reference, they will be
   * required often
   */

  for(AtomIndexType i = 0; i < _adjacencies.size(); i++) {
    for(const auto& j: _adjacencies.getAdjacencies(i)) {
      // avoid duplicate work by avoiding i > j cases
      if(i > j) continue;

      // enter the constraints for i, j
      const auto distance = Bond::calculateBondDistance(
        getElementType(i),
        getElementType(j),
        getBondType(i, j)
      );

      distanceBounds.upperBound(i, j) = (1 + oneTwoDelta) * distance;
      distanceBounds.lowerBound(i, j) = (1 - oneTwoDelta) * distance;
    }
  }

  return distanceBounds;
}

std::vector<DistanceConstraint> Molecule::_createConstraint(
  const std::vector<AtomIndexType>& chain
) const {
  auto& i = chain.front();
  auto& j = chain.back();

  switch(chain.size()) {
    case 2: {
      auto distance = Bond::calculateBondDistance(
        getElementType(i),
        getElementType(j),
        getBondType(i, j)
      );
      return {
        DistanceConstraint(
          i,
          j,
          0.99 * distance,
          1.01 * distance
        )
      };
    }
    case 3: {
      auto& intermediate = chain[1];
      // TODO continue here
      // no considerations of charge at all so far, not even VSEPR
      boost::optional<double> angle;
      // determine which symmetry applies at intermediate position, angle
      auto nLigands = getBondedAtomIndices(intermediate).size();
      if(nLigands == 2) {
        if(getElementType(intermediate) == Delib::ElementType::O) {
          angle = 104.5;
        } else {
          angle = Symmetry::angleFunction(Symmetry::Name::Linear)(0, 1);
        }
      } else if(nLigands == 3) {
        angle = Symmetry::angleFunction(Symmetry::Name::TrigonalPlanar)(0, 1);
      } else if(nLigands == 4) {
        angle = Symmetry::angleFunction(Symmetry::Name::Tetrahedral)(0, 1);
      } 

      if((bool) angle) {
        auto a = Bond::calculateBondDistance(
          getElementType(i),
          getElementType(intermediate),
          getBondType(i, intermediate)
        );
        auto b = Bond::calculateBondDistance(
          getElementType(intermediate),
          getElementType(j),
          getBondType(intermediate, j)
        );
        auto distance = CommonTrig::lawOfCosines(
          a,
          b,
          angle.value() * M_PI / 180.0
        );

        return {
          DistanceConstraint(
            i,
            j,
            0.95 * distance,
            1.05 * distance
          )
        };
      } else {
        return {};
      }
    }
    case 4: {
      // TODO implement
      return {};
    }
    default: {
      return {
        DistanceConstraint(
          i,
          j,
          (
            AtomInfo::vdwRadius(getElementType(i))
            + AtomInfo::vdwRadius(getElementType(j))
          ),
          100
        )
      };
    }
  }

}

std::vector<ChiralityConstraint> Molecule::getChiralityConstraints() const {
  std::vector<ChiralityConstraint> constraints;
  for(const auto& stereocenterPtr : _stereocenters) {
    for(const auto& distanceConstraint : stereocenterPtr -> chiralityConstraints() ) {
      constraints.push_back(distanceConstraint);
    }
  }
  return constraints;
}

} // eo namespace
