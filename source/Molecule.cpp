#include "boost/graph/biconnected_components.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/graph_utility.hpp"

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

  // Find CNStereocenters
  for(
    AtomIndexType candidateIndex = 0;
    candidateIndex < numAtoms();
    ++candidateIndex
  ) {
    if(_isCNStereocenterCandidate(candidateIndex)) {
      // Construct a Stereocenter here
      auto newStereocenter = std::make_shared<Stereocenters::CNStereocenter>(
        determineLocalGeometry(candidateIndex),
        candidateIndex,
        rankPriority(candidateIndex)
      );

      if(newStereocenter -> numAssignments() > 1) {
        stereocenterList.add(
          std::move(newStereocenter)
        );
      }
    }
  }

  return stereocenterList;
}


bool Molecule::_isValidIndex(const AtomIndexType& index) const {
  return index < numAtoms();
}

bool Molecule::_isCNStereocenterCandidate(
  const AtomIndexType& atomIndex,
  const TemperatureRegimeOption& temperatureRegime
) const {
  auto numAdjacencies = getNumAdjacencies(atomIndex);

  if(numAdjacencies < 3) {
    return false;
  }

  if(temperatureRegime == TemperatureRegimeOption::HighTemperature) {
    /* Skip any instances of nitrogen with exactly three adjacencies,
     * unless the central atom is part of a cycle of size 4 or smaller
     */
    if(
      getElementType(atomIndex) == Delib::ElementType::N 
      && numAdjacencies == 3
    ) {
      auto cycleData = getCycleData();

      // Find out if the nitrogen is in a cycle of size 4 or smaller
      bool isInCycleOfSize4OrSmaller = false;

      auto cycleIter = cycleData.getCyclesIteratorContaining(atomIndex);
      while(!cycleIter.atEnd()) {
        if(cycleIter.cycleSize() <= 4) {
          isInCycleOfSize4OrSmaller = true;
          break;
        }

        cycleIter.advance();
      }

      return isInCycleOfSize4OrSmaller;
    }
  }

  return true;
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

void Molecule::_pickyFitStereocenter(
  Stereocenters::CNStereocenter& stereocenter,
  const Symmetry::Name& expectedSymmetry,
  const Delib::PositionCollection& positions
) const {
  auto& centralAtom = stereocenter.getCentralAtomIndex();

  /* Seesaw and tetrahedral are surprisingly close in terms of angles, and
   * sometimes just slightly distorted tetrahedral centers can be recognized
   * as seesaws, even though it makes absolutely zero sense. So in case
   * the atom is a carbon, the expected geometry is tetrahedral and it has
   * four adjacencies, just exclude Seesaw from the list of symmetries being
   * fitted against.
   *
   * Calling
   * determineLocalGeometry is somewhat overkill here, but possibly more
   * future-proof.
   */
  if(
    getElementType(centralAtom) == Delib::ElementType::C
    && getNumAdjacencies(centralAtom) == 4
    && expectedSymmetry == Symmetry::Name::Tetrahedral
  ) {
    stereocenter.fit(
      positions,
      {Symmetry::Name::Seesaw}
    );
  } else {
    stereocenter.fit(positions);
  }
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

  // Ensure this is only called on non-terminal atoms
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
Molecule::Molecule() noexcept
  : Molecule(Delib::ElementType::H, Delib::ElementType::H, BondType::Single) {}

Molecule::Molecule(
  const Delib::ElementType& a,
  const Delib::ElementType& b,
  const BondType& bondType
) noexcept {
  // update _adjacencies
  _addAtom(a);
  _addAtom(b);
  // Although addBond is potentially-throwing, it never will
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
  if(!_isValidIndex(adjacentTo)) {
    throw std::out_of_range("Molecule::addAtom: Supplied atom index is invalid!");
  }

  const auto index = _addAtom(elementType);
  addBond(index, adjacentTo, bondType);
  return index;
}

void Molecule::addBond(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const BondType& bondType
) {
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::addBond: A supplied index is invalid!");
  }

  if(a == b) {
    throw std::logic_error("Molecule::addBond: Cannot add a bond between identical indices!");
  }

  auto edgeAddPair = boost::add_edge(a, b, _adjacencies);
  _adjacencies[edgeAddPair.first].bondType = bondType;
}

void Molecule::assignStereocenterAtAtom(
  const AtomIndexType& a,
  const boost::optional<unsigned>& assignment
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::assignStereocenterAtAtom: Supplied index is invalid!");
  }

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

void Molecule::refreshStereocenters() {
  _stereocenters = _detectStereocenters();
}

void Molecule::removeAtom(const AtomIndexType& a) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::removeAtom: Supplied index is invalid!");
  }

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
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::removeBond: Supplied index is invalid!");
  }

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

bool Molecule::setBondType(
  const AtomIndexType& a,
  const AtomIndexType& b,
  const BondType& bondType
) {
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::setBondType: A supplied index is invalid!");
  }

  auto edgePair = boost::edge(a, b, _adjacencies);
  if(edgePair.second) {
    _adjacencies[edgePair.first].bondType = bondType;
  } else {
    addBond(a, b, bondType);
  }

  return edgePair.second;
}

void Molecule::setElementType(
  const AtomIndexType& a,
  const Delib::ElementType& elementType
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::setElementType: This index is invalid!");
  }

  _adjacencies[a].elementType = elementType;
}

void Molecule::setGeometryAtAtom(
  const AtomIndexType& a,
  const Symmetry::Name& symmetryName
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::setGeometryAtAtom: Supplied atom index is invalid");
  }

  if(_stereocenters.involving(a)) {
    if(_stereocenters.at(a)->type() == Stereocenters::Type::CNStereocenter) {
      auto CNSPointer = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
        _stereocenters.at(a)
      );

      if(
        Symmetry::size(CNSPointer->getSymmetry()) 
        == Symmetry::size(symmetryName)
      ) {
        CNSPointer->setSymmetry(symmetryName);
      } else {
        throw std::logic_error(
          "Molecule::setGeometryAtAtom: The size of the supplied symmetry is "
          "not the same as that of the existing stereocenter's current symmetry!"
        );
      }
    } else {
      throw std::logic_error(
        "Molecule::setGeometryAtAtom: There is an E/Z stereocenter at the "
        "supplied position, its local symmetry cannot be changed"
      );
    }
  } else {
    const Symmetry::Name expectedSymmetry = determineLocalGeometry(a);

    if(Symmetry::size(expectedSymmetry) != Symmetry::size(symmetryName)) {
      throw std::logic_error(
        "Molecule::setGeometryAtAtom: The size of the supplied symmetry is not "
        " the same as that of the expected symmetry!"
      );
    }

    /* In case the expected symmetry is the same as the supplied symmetry, this
     * is a no-op, since it adds no information to the Molecule (this
     * stereocenter would otherwise be created during DG initialization)
     */
    if(expectedSymmetry != symmetryName) {
      // Add the stereocenter irrespective of how many assignments it has
      auto newStereocenterPtr = std::make_shared<Stereocenters::CNStereocenter>(
        symmetryName,
        a,
        rankPriority(a)
      );

      // Default-assign stereocenters with only one assignment
      if(newStereocenterPtr->numAssignments() == 1) {
        newStereocenterPtr->assign(0u);
      }

      _stereocenters.add(
        std::move(newStereocenterPtr)
      );
    }
  }
}

/* Information */
Symmetry::Name Molecule::determineLocalGeometry(
  const AtomIndexType& index
) const {
  if(!_isValidIndex(index)) {
    throw std::out_of_range("Molecule::determineLocalGeometry: Supplied index is invalid!");
  }

  if(getNumAdjacencies(index) <= 1) {
    throw std::logic_error(
      "Molecule::determineLocalGeometry: No geometries exist for terminal or "
      " disconnected atoms"
    );
  }

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
    edges.push_back(
      ExplicitEdge({
        {boost::source(edgeIndex, _adjacencies), boost::target(edgeIndex, _adjacencies)},
        _adjacencies[edgeIndex].bondType
      })
    );
  }

  return edges;
}

Delib::ElementType Molecule::getElementType(const AtomIndexType& index) const {
  if(!_isValidIndex(index)) {
    throw std::out_of_range("Molecule::getElementType: Supplied index is invalid");
  }

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

  for(const auto& edgeIndex : _getEZStereocenterCandidates()) {
    auto source = boost::source(edgeIndex, _adjacencies),
         target = boost::target(edgeIndex, _adjacencies);

    // Construct a Stereocenter here
    auto newStereocenter = std::make_shared<Stereocenters::EZStereocenter>(
      source,
      rankPriority(source, {target}, positions),
      target,
      rankPriority(target, {source}, positions)
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
    // Skip unsuitable atoms and ones that already have a stereocenter
    if(
      !_isCNStereocenterCandidate(candidateIndex)
      || stereocenters.involving(candidateIndex)
    ) {
      continue;
    }

    const Symmetry::Name expectedGeometry = determineLocalGeometry(candidateIndex);

    // Construct it
    auto stereocenterPtr = std::make_shared<Stereocenters::CNStereocenter>(
      expectedGeometry,
      candidateIndex,
      rankPriority(candidateIndex, {}, positions)
    );

    _pickyFitStereocenter(
      *stereocenterPtr,
      expectedGeometry,
      positions
    );

    /* Add the CNStereocenter to the list, unless if it has one assignment only
     * and the symmetry is the same as expected
     */
    if(
      !(
        stereocenterPtr -> numAssignments() == 1
        && expectedGeometry == stereocenterPtr -> getSymmetry()
      )
    ) {
      stereocenters.add(
        std::move(stereocenterPtr)
      );
    }
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
  // A molecule is by definition at least two atoms!
  if(numAtoms() == 2) {
    return false;
  }

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
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::iterateAdjacencies: Supplied index is invalid!");
  }

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
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::iterateEdges: Supplied index is invalid!");
  }

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
  const std::set<AtomIndexType>& excludeAdjacent,
  const boost::optional<Delib::PositionCollection>& positionsOption
) const {
  RankingInformation rankingResult;

  // Rank the substituents
  auto expandedTree = RankingTree(
    *this,
    a,
    excludeAdjacent,
    RankingTree::ExpansionOption::Optimized,
    positionsOption
  );

  rankingResult.sortedSubstituents = expandedTree.getRanked();

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
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::operator[]: Supplied index is invalid!");
  }
  
  return RangeForTemporary<GraphType::adjacency_iterator>(
    boost::adjacent_vertices(a, _adjacencies)
  );
}

bool Molecule::operator == (const Molecule& other) const {
  const unsigned thisNumAtoms = numAtoms();

  if(thisNumAtoms != other.numAtoms()) {
    return false;
  }

  // Where the corresponding index from the other graph is stored
  std::vector<AtomIndexType> indexMap(numAtoms());

  bool isomorphic = boost::isomorphism(
    _adjacencies,
    other._adjacencies,
    boost::isomorphism_map(
      boost::make_safe_iterator_property_map(
        indexMap.begin(),
        thisNumAtoms,
        boost::get(boost::vertex_index, _adjacencies)
      )
    )
  );

  if(!isomorphic) {
    return false;
  }

  // Check that all element types of the isomorphism mapping match
  for(AtomIndexType thisIndex = 0; thisIndex < thisNumAtoms; ++thisIndex) {
    if(
      _adjacencies[thisIndex].elementType 
        != other._adjacencies[indexMap.at(thisIndex)].elementType
    ) {
      return false;
    }
  }

  // Check that all bond types are identical
  for(const auto& edgeIndex : iterateEdges()) {
    // Fetch the corresponding edge from the other graph
    auto otherEdgePair = boost::edge(
      indexMap.at(boost::source(edgeIndex, _adjacencies)),
      indexMap.at(boost::target(edgeIndex, _adjacencies)),
      other._adjacencies
    );

    // This edge MUST be found, the isomorphism holds
    assert(otherEdgePair.second);

    if(
      _adjacencies[edgeIndex].bondType
      != other._adjacencies[otherEdgePair.first].bondType
    ) {
      return false;
    }
  }

  // Before doing a full equivalence check, peek at the sizes
  if(_stereocenters.size() != other._stereocenters.size()) {
    return false;
  }

  // Check equivalence of the StereocenterLists
  for(const auto& stereocenterPtr : _stereocenters) {
    if(stereocenterPtr->type() == Stereocenters::Type::CNStereocenter) {
      const auto otherCentralAtom = indexMap.at(
        stereocenterPtr->involvedAtoms().front()
      );

      // Does other not have a stereocenter there?
      if(!other._stereocenters.involving(otherCentralAtom)) {
        return false;
      }

      // Ensure the type at other is a CNStereocenter too
      if(
        other._stereocenters.at(otherCentralAtom)->type() 
          != Stereocenters::Type::CNStereocenter
      ) {
        return false;
      }

      auto thisCNSPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(stereocenterPtr);
      auto otherCNSPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
        other._stereocenters.at(otherCentralAtom)
      );

      // Are they equal in an abstract sense, not object-representation-wise?
      if(
        thisCNSPtr->getSymmetry() != otherCNSPtr->getSymmetry()
        || thisCNSPtr->numAssignments() != otherCNSPtr->numAssignments()
        || thisCNSPtr->assigned() != otherCNSPtr->assigned()
      ) {
        return false;
      }
    } else {
      const auto otherCentralAtoms = TemplateMagic::map(
        stereocenterPtr->involvedAtoms(),
        [&indexMap](const auto& thisIndex) -> AtomIndexType {
          return indexMap.at(thisIndex);
        }
      );

      // Ensure there is a stereocenter on both
      if(
        !other._stereocenters.involving(otherCentralAtoms.front())
        || !other._stereocenters.involving(otherCentralAtoms.back())
      ) {
        return false;
      }

      // Address-compare that the stereocenters on other are identical for both
      if(
        other._stereocenters.at(otherCentralAtoms.front()) 
          != other._stereocenters.at(otherCentralAtoms.back())
      ) {
        return false;
      }

      // Ensure that it's also an EZStereocenter
      if(
        other._stereocenters.at(otherCentralAtoms.front())->type() 
          != Stereocenters::Type::EZStereocenter
      ) {
        return false;
      }

      // Shortcut name
      const auto& otherPtr = other._stereocenters.at(otherCentralAtoms.front());

      // Abstract-compare both
      if(
        stereocenterPtr->numAssignments() != otherPtr->numAssignments()
        || stereocenterPtr->assigned() != otherPtr->assigned()
      ) {
        return false;
      }
    }
  }

  return true;
}

bool Molecule::operator != (const Molecule& other) const {
  return !(*this == other);
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
