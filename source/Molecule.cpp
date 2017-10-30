#include "boost/graph/biconnected_components.hpp"
#include "boost/graph/graphviz.hpp"

#include "template_magic/Containers.h"
#include "template_magic/Numeric.h"

#include "CNStereocenter.h"
#include "CommonTrig.h"
#include "EZStereocenter.h"
#include "GraphAlgorithms.h"
#include "Log.h"
#include "Molecule.h"
#include "MolGraphWriter.h"
#include "RankingTree.h"

namespace MoleculeManip {

/* Molecule implementation ---------------------------------------------------*/
/* Private members */

AtomIndexType Molecule::_addAtom(const Delib::ElementType& elementType) {
  auto vertex = boost::add_vertex(_adjacencies);
  _adjacencies[vertex].elementType = elementType;
  return vertex;
}

StereocenterList Molecule::_detectStereocenters() const {
  StereocenterList stereocenterList;

  // Find CNStereocenters
  for(const auto& candidateIndex : _getCNStereocenterCandidates()) {
    // Construct a Stereocenter here
    auto newStereocenter = std::make_shared<
      Stereocenters::CNStereocenter
    >(
      determineLocalGeometry(candidateIndex),
      candidateIndex,
      rankPriority(candidateIndex)
    );

    /*std::cout << "Trial stereocenter: " << newStereocenter -> info() << std::endl;
    std::cout << "Ranked adjacent indices (low to high): vec{";
    for(const auto& adjacency: rankResultPair.first) {
      std::cout << adjacency;
      if(adjacency != rankResultPair.first.back()) std::cout << ", ";
    }
    std::cout << "}" << std::endl;

    std::cout << "Equal pairs: vec{";
    for(const auto& indexPair : rankResultPair.second) {
      std::cout << "(" << indexPair.first << ", " << indexPair.second << ")";
    }
    std::cout << "}" << std::endl;*/

    if(newStereocenter -> numAssignments() > 1) {
      stereocenterList.add(
        std::move(newStereocenter)
      );
    }
  }

  /* TODO
   * - Will need refinement to not instantiate EZStereocenters in small cycles
   *   (up to a preset size, maybe around 8 or so?)
   */
  // Find EZStereocenters
  for(const auto& edgeIndex : _getEZStereocenterCandidates()) {
    auto source = boost::source(edgeIndex, _adjacencies),
         target = boost::target(edgeIndex, _adjacencies);

    // Construct a Stereocenter here
    auto newStereocenter = std::make_shared<
      Stereocenters::EZStereocenter
    >(
      source,
      rankPriority(source, {target}),
      target,
      rankPriority(target, {source})
    );

    if(newStereocenter -> numAssignments() == 2) {
      stereocenterList.add(
        std::move(newStereocenter)
      );
    }
  }

  return stereocenterList;
}


bool Molecule::_isValidIndex(const AtomIndexType& index) const {
  return index < numAtoms();
}

std::vector<AtomIndexType> Molecule::_getCNStereocenterCandidates() const {
  std::vector<AtomIndexType> candidates;

  for(AtomIndexType i = 0; i < numAtoms(); i++) {
    if(
      /* TODO this is no longer a valid way of checking how many ligands there are
       * -> eta bonds exist!
       */
      getNumAdjacencies(i) >= 3 
    ) {
      candidates.push_back(i);
    }
  }

  return candidates;
}

std::vector<EdgeIndexType> Molecule::_getEZStereocenterCandidates() const {
  auto numNonEtaAdjacencies = [&](const AtomIndexType& a) -> unsigned {
    unsigned nonEta = 0;

    for(const auto& edgeIndex : iterateEdges(a)) {
      if(_adjacencies[edgeIndex].bondType != BondType::Eta) {
        nonEta += 1;
      }
    }

    return nonEta;
  };

  std::vector<EdgeIndexType> candidates;

  for(
    const auto& edgeIndex : 
    RangeForTemporary<GraphType::edge_iterator>(
      boost::edges(_adjacencies)
    )
  ) {
    auto source = boost::source(edgeIndex, _adjacencies),
         target = boost::target(edgeIndex, _adjacencies);

    auto sourceAdjacencies = numNonEtaAdjacencies(source),
         targetAdjacencies = numNonEtaAdjacencies(target);

    if(
      _adjacencies[edgeIndex].bondType == BondType::Double
      && 2 <= sourceAdjacencies
      && sourceAdjacencies <= 3
      && 2 <= targetAdjacencies 
      && targetAdjacencies <= 3
    ) {
      candidates.push_back(edgeIndex);
    }
  }

  return candidates;
}

std::vector<LocalGeometry::LigandType> Molecule::_reduceToLigandTypes(
  const AtomIndexType& index
) const {
  /* TODO 
   * - No L, X determination. Although, will L, X even be needed for metals?
   *   Maybe only for OZ and NVE determination...
   */
  /* VSEPR formulation is that geometry is a function of 
   * - localized charge of central atom
   * - atom type of central atom, neighbors
   * - bond types to neighbors
   */

  // call this only on non-terminal atoms
  assert(getNumAdjacencies(index) > 1);

  // first basic stuff for VSEPR, later L and X for transition metals
  // geometry inference does not care if the substituents are somehow 
  // connected (unless in later models the entire structure is considered)
  std::vector<LocalGeometry::LigandType> ligands;

  for(const auto& adjacentIndex: iterateAdjacencies(index)) {
    ligands.push_back(
      LocalGeometry::LigandType {
        0, 0, {
          {  // L and X are 0 since only VSEPR is considered for now
            getElementType(adjacentIndex),
            getBondType(index, adjacentIndex).value()
          }
        }
      }
    );
  }

  return ligands;
}

void Molecule::_updateStereocenterList() {
  /* Two cases: If the StereocenterList is empty, we can just use detect to
   * find stereocenters in the Molecule, if there are any new ones.
   *
   * In the other case, we have to recheck every existing stereocenter. If
   * ranking was affected and the stereocenter has a set assignment, we need to
   * find the assignment that the previous ranking represented spatially in the
   * new set of assignments and assign the stereocenter to that.
   */
  if(_stereocenters.empty()) {
    _stereocenters = _detectStereocenters();
  } else {
    // TODO this is difficult!
  }
}


/* Public members */
/* Constructors */
Molecule::Molecule(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) {
  // update _adjacencies
  _addAtom(a);
  _addAtom(b);
  addBond(0, 1, bondType);
}

Molecule::Molecule(
  const Delib::ElementTypeCollection& elements,
  const Edges& edges
) {
  for(const auto& element: elements) {
    _addAtom(element);
  }

  for(const auto& edge: edges) {
    addBond(edge.first.first, edge.first.second, edge.second);
  }
}

Molecule::Molecule(const GraphType& graph) 
: _adjacencies(graph), 
  _stereocenters(_detectStereocenters())
{}

Molecule::Molecule(
  const GraphType& graph,
  const Delib::PositionCollection& positions
) : _adjacencies(graph),
    _stereocenters(inferStereocentersFromPositions(positions))
{}

/* Modifiers */
AtomIndexType Molecule::addAtom(
  const Delib::ElementType& elementType,
  const AtomIndexType& adjacentTo,
  const BondType& bondType
) {
  const auto index = _addAtom(elementType);
  addBond(index, adjacentTo, bondType);
  return index;
}

void Molecule::addBond(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const BondType& bondType
) {
  assert(_isValidIndex(a) && _isValidIndex(b) && a != b);
  auto edgeAddPair = boost::add_edge(a, b, _adjacencies);
  _adjacencies[edgeAddPair.first].bondType = bondType;
}

void Molecule::assignStereocenterAtAtom(
  const AtomIndexType& a,
  const boost::optional<unsigned>& assignment
) {
  assert(_isValidIndex(a));
  if(_stereocenters.involving(a)) {
    auto stereocenterPtr = _stereocenters.at(a);

    if(assignment < stereocenterPtr -> numAssignments()) {
      stereocenterPtr -> assign(assignment);

      // TODO keep valid state trigger, changing an assignment can alter ranking
    } else {
      throw std::logic_error("assignStereocenterAtAtom: Invalid assignment index!");
    }
  } else {
    throw std::logic_error("assignStereocenterAtAtom: No stereocenter at this index!");
  }
}

void Molecule::changeElementType(
  const AtomIndexType& a,
  const Delib::ElementType& elementType
) {
  if(!_isValidIndex(a)) {
    throw std::logic_error("This index is invalid!");
  }

  _adjacencies[a].elementType = elementType;
}

void Molecule::refreshStereocenters() {
  _stereocenters = _detectStereocenters();
}

void Molecule::removeAtom(const AtomIndexType& a) {
  if(!isSafeToRemoveAtom(a)) {
    throw std::logic_error("Removing this atom disconnects the graph!");
  }

  // Remove all edges to and from this vertex
  boost::clear_vertex(a, _adjacencies);

  // Remove the vertex itself
  boost::remove_vertex(a, _adjacencies);
}

void Molecule::removeBond(
  const AtomIndexType& a,
  const AtomIndexType& b
) {
  if(!isSafeToRemoveBond(a, b)) {
    throw std::logic_error("Removing this bond disconnects the graph!");
  }

  // Find edge
  auto edgePair = boost::edge(a, b, _adjacencies);
  if(edgePair.second) {
    boost::remove_edge(edgePair.first, _adjacencies);
  } else {
    throw std::logic_error("That bond does not exist!");
  }
}

/* Information */
Symmetry::Name Molecule::determineLocalGeometry(
  const AtomIndexType& index
) const {
  assert(getNumAdjacencies(index) > 1); 

  auto ligandsVector = _reduceToLigandTypes(index);

  // TODO this below is invalid for metals!
  unsigned nSites = getNumAdjacencies(index);
  int formalCharge = 0;

  if(AtomInfo::isMainGroupElement(getElementType(index))) {
    return LocalGeometry::VSEPR::determineGeometry(
      getElementType(index),
      nSites,
      ligandsVector,
      formalCharge
    );
  } 

  // Pick the first Symmetry of fitting size
  auto findIter = std::find_if(
    Symmetry::allNames.begin(),
    Symmetry::allNames.end(),
    [&nSites](const auto& symmetryName) -> bool {
      return Symmetry::size(symmetryName) == nSites;
    }
  );

  if(findIter == Symmetry::allNames.end()) {
    throw std::logic_error(
      "Could not find a suitable local geometry!"
    );
  }

  return *findIter;
}

std::string Molecule::dumpGraphviz() const {
  MolGraphWriter propertyWriter(&_adjacencies);

  std::stringstream graphvizStream;

  boost::write_graphviz(
    graphvizStream,
    _adjacencies,
    propertyWriter,
    propertyWriter,
    propertyWriter
  );

  return graphvizStream.str();
}

std::vector<AtomIndexType> Molecule::getAdjacencies(
  const AtomIndexType& a
) const {
  std::vector<AtomIndexType> copy;

  // C++17 auto [begin, end] = ...
  GraphType::adjacency_iterator begin, end;
  std::tie(begin, end) = boost::adjacent_vertices(a, _adjacencies);
  std::copy(
    begin,
    end,
    std::back_inserter(copy)
  );

  return copy;
}

boost::optional<BondType> Molecule::getBondType(
  const AtomIndexType& a,
  const AtomIndexType& b
) const {
  auto edgePair = boost::edge(a, b, _adjacencies);

  if(edgePair.second) {
    return _adjacencies[edgePair.first].bondType;
  } 

  // fallback
  return boost::none;
}

CycleData Molecule::getCycleData() const {
  return CycleData(_adjacencies);
}

// Creates a copy of the contained data suitable for the Edges class
std::vector<Molecule::ExplicitEdge> Molecule::getEdges() const {
  std::vector<ExplicitEdge> edges;

  for(
    const auto& edgeIndex : 
    RangeForTemporary<GraphType::edge_iterator>(boost::edges(_adjacencies))
  ) {
    edges.push_back(ExplicitEdge({
      {boost::source(edgeIndex, _adjacencies), boost::target(edgeIndex, _adjacencies)},
      _adjacencies[edgeIndex].bondType
    }));
  }

  return edges;
}

Delib::ElementType Molecule::getElementType(const AtomIndexType& index) const {
  assert(_isValidIndex(index));
  return _adjacencies[index].elementType;
}

const GraphType& Molecule::getGraph() const {
  return _adjacencies;
}

const StereocenterList& Molecule::getStereocenterList() const {
  return _stereocenters;
}

unsigned Molecule::getNumAdjacencies(
  const AtomIndexType& a
) const {
  return boost::out_degree(a, _adjacencies);
}

StereocenterList Molecule::inferStereocentersFromPositions(
  const Delib::PositionCollection& positions
) const {
  StereocenterList stereocenters;

  for(
    const auto& edgeIndex : 
    _getEZStereocenterCandidates()
  ) {
    auto source = boost::source(edgeIndex, _adjacencies),
         target = boost::target(edgeIndex, _adjacencies);

    // Construct a Stereocenter here
    auto newStereocenter = std::make_shared<
      Stereocenters::EZStereocenter
    >(
      source,
      rankPriority(source, {target}),
      target,
      rankPriority(target, {source})
    );

    newStereocenter -> fit(positions);

    if(newStereocenter -> numAssignments() == 2) {
      stereocenters.add(
        std::move(newStereocenter)
      );
    }
  }

  /* Add a CNStereocenter everywhere where the symmetry yielding the best fit is 
   * not the one that Molecule's determineLocalGeometry gets and where we
   * can fully determine a Stereocenter's assignment from the positions
   */
  for(unsigned candidateIndex = 0; candidateIndex < numAtoms(); candidateIndex++) {
    // Skip terminal atoms and ones that already have a stereocenter
    if(
      getNumAdjacencies(candidateIndex) <= 1
      || stereocenters.involving(candidateIndex)
    ) {
      continue;
    }

    // Construct it
    std::shared_ptr<
      Stereocenters::CNStereocenter
    > stereocenterPtr = std::make_shared<
      Stereocenters::CNStereocenter
    >(
      determineLocalGeometry(candidateIndex),
      candidateIndex,
      rankPriority(candidateIndex)
    );

    stereocenterPtr -> fit(positions);

    /* TODO In case the CNStereocenter has one assignment only and the symmetry
     * is the same as by determined, do not add it to the list of stereocenters
     */
    stereocenters.add(
      std::move(stereocenterPtr)
    );
  }

  // TODO EZStereocenters
  /* NOTES
   * - CNStereocenter detection may have generated trigonal planar
   *   stereocenters on the endpoints of the double bond edge -> remove if an
   *   EZStereocenter is instantiated there instead
   * - Issues: double bond and CNStereocenter can clash -> Grubbs
   *
   * STEPS
   * - Calculate dihedral angle of high-priority pair from 3D
   *   -> Select E/Z within tolerance of 0° / 180° endpoints
   *   -> Throw outside of those tolerances
   */

  return stereocenters;
}


bool Molecule::isAdjacent(
  const AtomIndexType& a,
  const AtomIndexType& b
) const {
  GraphType::adjacency_iterator begin, end;
  std::tie(begin, end) = boost::adjacent_vertices(a, _adjacencies);
  return std::find(
    begin,
    end,
    b
  ) != end;
}

bool Molecule::isSafeToRemoveAtom(const AtomIndexType& a) const {
  auto removalSafetyData = GraphAlgorithms::getRemovalSafetyData(
    getGraph()
  );

  return removalSafetyData.articulationVertices.count(a) == 0;
}

bool Molecule::isSafeToRemoveBond(
  const AtomIndexType& a,
  const AtomIndexType& b
) const {
  EdgeIndexType edgeIndex;
  bool foundEdge;

  std::tie(edgeIndex, foundEdge) = boost::edge(a, b, _adjacencies);

  if(foundEdge) {
    auto removalSafetyData = GraphAlgorithms::getRemovalSafetyData(
      getGraph()
    );

    return removalSafetyData.bridges.count(edgeIndex) == 0;
  }

  // Bond does not exist -> it is unsafe to remove
  return false;
}


/*! Returns a range-for temporary object allowing c++11 style for loop 
 * iteration through an atom's adjacencies
 */
RangeForTemporary<GraphType::adjacency_iterator> Molecule::iterateAdjacencies(
  const AtomIndexType& a
) const {
  assert(_isValidIndex(a));

  return RangeForTemporary<GraphType::adjacency_iterator>(
    boost::adjacent_vertices(a, _adjacencies)
  );
}

RangeForTemporary<GraphType::edge_iterator> Molecule::iterateEdges() const {
  return RangeForTemporary<GraphType::edge_iterator>(boost::edges(_adjacencies));
}

RangeForTemporary<GraphType::out_edge_iterator> Molecule::iterateEdges(
  const AtomIndexType& a
) const {
  assert(_isValidIndex(a));

  return RangeForTemporary<GraphType::out_edge_iterator>(
    boost::out_edges(a, _adjacencies)
  );
}

unsigned Molecule::numAtoms() const {
  return boost::num_vertices(_adjacencies);
}

unsigned Molecule::numBonds() const {
  return boost::num_edges(_adjacencies);
}

RankingInformation Molecule::rankPriority(
  const AtomIndexType& a,
  const std::set<AtomIndexType>& excludeAdjacent 
) const {
  RankingInformation rankingResult;

  // Rank the substituents
  auto expandedTree = RankingTree(*this, a);
  rankingResult.sortedSubstituents = expandedTree.rank(excludeAdjacent);

  auto adjacentIndices = getAdjacencies(a);
  std::set<AtomIndexType> activeIndices;
  for(const auto& adjacentIndex : adjacentIndices) {
    if(excludeAdjacent.count(adjacentIndex) == 0) {
      activeIndices.insert(adjacentIndex);
    }
  }


  // Find links between them
  rankingResult.linkedPairs = GraphAlgorithms::findSubstituentLinks(
    _adjacencies,
    a,
    activeIndices
  );

  return rankingResult;
}

/* Operators */
RangeForTemporary<GraphType::adjacency_iterator> Molecule::operator [] (
  const AtomIndexType& a
) const {
  assert(_isValidIndex(a));
  
  return RangeForTemporary<GraphType::adjacency_iterator>(
    boost::adjacent_vertices(a, _adjacencies)
  );
}

} // namespace MoleculeManip

std::ostream& operator << (
  std::ostream& os,
  const MoleculeManip::Molecule& molecule
) {
  if(!molecule.getStereocenterList().empty()) {
    os << "Stereocenter information:\n";

    for(const auto& stereocenterPtr: molecule.getStereocenterList()) {
      os << stereocenterPtr -> info() << std::endl;
    }
  }

  return os;
}
