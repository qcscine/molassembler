#include "molassembler/Molecule/MoleculeImpl.h"

#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/graph_utility.hpp"
#include "chemical_symmetries/ConstexprProperties.h"
#include "Delib/Constants.h"
#include "Delib/ElementTypeCollection.h"

#include "molassembler/Cycles.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Modeling/LocalGeometryModel.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/Molecule/MolGraphWriter.h"
#include "molassembler/Molecule/RankingTree.h"
#include "molassembler/Options.h"

namespace molassembler {

StereocenterList Molecule::Impl::_detectStereocenters() const {
  StereocenterList stereocenterList;

  Cycles cycleData = graph().cycles();
  // Find AtomStereocenters
  for(
    AtomIndex candidateIndex = 0;
    candidateIndex < graph().N();
    ++candidateIndex
  ) {
    RankingInformation localRanking = rankPriority(candidateIndex);

    // Skip terminal atoms
    if(localRanking.ligands.size() <= 1) {
      continue;
    }

    // Construct a Stereocenter here
    auto newStereocenter = AtomStereocenter {
      _adjacencies,
      determineLocalGeometry(candidateIndex, localRanking),
      candidateIndex,
      localRanking
    };

    if(newStereocenter.numAssignments() == 1) {
      newStereocenter.assign(0);
    }

    if(
      !disregardStereocenter(
        newStereocenter,
        graph().elementType(candidateIndex),
        cycleData,
        Options::temperatureRegime
      )
    ) {
      stereocenterList.add(
        candidateIndex,
        std::move(newStereocenter)
      );
    }
  }

  const InnerGraph& inner = graph().inner();

  // Find BondStereocenters
  for(
    const InnerGraph::Edge& edgeIndex :
    boost::make_iterator_range(inner.edges())
  ) {
    // Skip edges that cannot be candidates
    if(!_isGraphBasedBondStereocenterCandidate(edgeIndex)) {
      continue;
    }

    AtomIndex source = inner.source(edgeIndex),
              target = inner.target(edgeIndex);

    auto sourceAtomStereocenterOption = stereocenterList.option(source);
    auto targetAtomStereocenterOption = stereocenterList.option(target);

    // There need to be assigned stereocenters on both vertices
    if(
      !sourceAtomStereocenterOption
      || !targetAtomStereocenterOption
      || sourceAtomStereocenterOption->assigned() == boost::none
      || targetAtomStereocenterOption->assigned() == boost::none
    ) {
      continue;
    }

    const auto outerEdgeIndex = BondIndex {source, target};

    // Construct a Stereocenter here
    auto newStereocenter = BondStereocenter {
      *sourceAtomStereocenterOption,
      *targetAtomStereocenterOption,
      outerEdgeIndex
    };

    if(newStereocenter.numAssignments() > 1) {
      stereocenterList.add(
        outerEdgeIndex,
        std::move(newStereocenter)
      );
    }
  }

  return stereocenterList;
}

void Molecule::Impl::_ensureModelInvariants() const {
  if(graph().inner().connectedComponents() > 1) {
    throw std::logic_error("Molecules must be a single connected component. The supplied graph has multiple");
  }

  if(graph().N() < 2) {
    throw std::logic_error("Molecules must consist of at least two atoms!");
  }
}

bool Molecule::Impl::_isValidIndex(const AtomIndex index) const {
  return index < graph().N();
}

bool Molecule::Impl::_isGraphBasedBondStereocenterCandidate(
  const InnerGraph::Edge& e
) const {
  const BondType bondType = graph().inner().bondType(e);
  return (
    bondType == BondType::Double
    || bondType == BondType::Triple
    || bondType == BondType::Quadruple
    || bondType == BondType::Quintuple
    || bondType == BondType::Sextuple
  );
}

void Molecule::Impl::_propagateGraphChange() {
  /* Two cases: If the StereocenterList is empty, we can just use detect to
   * find any new stereocenters in the Molecule.
   */

  if(_stereocenters.empty()) {
    _stereocenters = _detectStereocenters();
    return;
  }

  /* In the other case, we have to recheck absolutely everywhere. If ranking
   * was affected and the stereocenter has a set assignment, we need to find
   * the assignment that the previous ranking represented spatially in the new
   * set of assignments and assign the stereocenter to that.
   */

  InnerGraph& inner = _adjacencies.inner();

  GraphAlgorithms::findAndSetEtaBonds(inner);

  /* TODO
   * - Need state propagation for BondStereocenters, anything else is madness
   */

  Cycles cycleData = graph().cycles();

  for(
    const InnerGraph::Vertex vertex :
    boost::make_iterator_range(inner.vertices())
  ) {
    auto stereocenterOption = _stereocenters.option(vertex);
    RankingInformation localRanking = rankPriority(vertex);

    if(stereocenterOption) {
      // The atom has become terminal
      if(localRanking.ligands.size() <= 1) {
        _stereocenters.remove(vertex);
        continue;
      }

      // Has the ranking changed?
      if(localRanking == stereocenterOption->getRanking()) {
        continue;
      }

      // Propagate the stereocenter state to the new ranking
      stereocenterOption -> propagateGraphChange(
        _adjacencies,
        localRanking
      );

      /* If the modified stereocenter has only one assignment and is unassigned
       * due to the graph change, default-assign it
       */
      if(
        stereocenterOption -> numStereopermutations() == 1
        && stereocenterOption -> numAssignments() == 1
      ) {
        stereocenterOption -> assign(0);
      }

      // If the change makes the stereocenter undesirable, remove it
      if(
        disregardStereocenter(
          *stereocenterOption,
          inner.elementType(vertex),
          cycleData,
          Options::temperatureRegime
        )
      ) {
        _stereocenters.remove(vertex);
      }

      /* If there are any BondStereocenters on adjacent edges, remove them
       * (since the Ranking has changed on a constituting AtomStereocenter)
       */
      for(
        const BondIndex& bond :
        boost::make_iterator_range(_adjacencies.bonds(vertex))
      ) {
        _stereocenters.try_remove(bond);
      }
    } else {
      // Skip terminal atoms
      if(localRanking.ligands.size() <= 1) {
        continue;
      }

      auto newStereocenter = AtomStereocenter {
        _adjacencies,
        determineLocalGeometry(vertex, localRanking),
        vertex,
        localRanking
      };

      // Default-assign trivial stereocenters
      if(
        newStereocenter.numStereopermutations() == 1
        && newStereocenter.numAssignments() == 1
      ) {
        newStereocenter.assign(0);
      }

      if(
        !disregardStereocenter(
          newStereocenter,
          inner.elementType(vertex),
          cycleData,
          Options::temperatureRegime
        )
      ) {
        _stereocenters.add(vertex, newStereocenter);
      }
    }
  }

  /* Any BondStereocenters whose constituing AtomStereocenter's rankings have
   * changed have been removed. So, I think now we just have to check if there
   * are any new ones we could have.
   */

  for(
    const InnerGraph::Edge& edge :
    boost::make_iterator_range(inner.edges())
  ) {
    // Check if the edge could be a stereocenter on graph based grounds
    if(!_isGraphBasedBondStereocenterCandidate(edge)) {
      continue;
    }

    // Fetch some indices
    AtomIndex source = inner.source(edge),
              target = inner.target(edge);

    BondIndex bond {source, target};

    // Ensure there isn't already a stereocenter on this edge
    if(auto stereocenterOption = _stereocenters.option(bond)) {
      continue;
    }

    // From here just like _detectStereocenters()

    auto sourceAtomStereocenterOption = _stereocenters.option(source);
    auto targetAtomStereocenterOption = _stereocenters.option(target);

    // There need to be assigned stereocenters on both vertices
    if(
      !sourceAtomStereocenterOption
      || !targetAtomStereocenterOption
      || sourceAtomStereocenterOption->assigned() == boost::none
      || targetAtomStereocenterOption->assigned() == boost::none
    ) {
      continue;
    }

    // Construct a Stereocenter here
    auto newStereocenter = BondStereocenter {
      *sourceAtomStereocenterOption,
      *targetAtomStereocenterOption,
      bond
    };

    if(newStereocenter.numAssignments() > 1) {
      _stereocenters.add(
        bond,
        std::move(newStereocenter)
      );
    }
  }
}

/* Public members */
/* Constructors */
Molecule::Impl::Impl() noexcept
  : Impl(Delib::ElementType::H, Delib::ElementType::H, BondType::Single) {}

Molecule::Impl::Impl(
  const Delib::ElementType a,
  const Delib::ElementType b,
  const BondType bondType
) noexcept {
  // update _adjacencies
  InnerGraph& inner = _adjacencies.inner();
  InnerGraph::Vertex i = inner.addVertex(a);
  InnerGraph::Vertex j = inner.addVertex(b);
  inner.addEdge(i, j, bondType);
}

Molecule::Impl::Impl(OuterGraph graph)
  : _adjacencies(std::move(graph))
{
  // Initialization
  GraphAlgorithms::findAndSetEtaBonds(_adjacencies.inner());
  _stereocenters = _detectStereocenters();
  _ensureModelInvariants();
}

Molecule::Impl::Impl(
  OuterGraph graph,
  const AngstromWrapper& positions,
  const boost::optional<
    std::vector<BondIndex>
  >& bondStereocenterCandidatesOptional
) : _adjacencies(std::move(graph))
{
  GraphAlgorithms::findAndSetEtaBonds(_adjacencies.inner());
  _stereocenters = inferStereocentersFromPositions(
    positions,
    bondStereocenterCandidatesOptional
  );
  _ensureModelInvariants();
}

Molecule::Impl::Impl(
  OuterGraph graph,
  StereocenterList stereocenters
) : _adjacencies(std::move(graph)),
    _stereocenters(std::move(stereocenters))
{
  // Initialization
  _ensureModelInvariants();
}

/* Modifiers */
AtomIndex Molecule::Impl::addAtom(
  const Delib::ElementType elementType,
  const AtomIndex adjacentTo,
  const BondType bondType
) {
  if(!_isValidIndex(adjacentTo)) {
    throw std::out_of_range("Molecule::addAtom: Supplied atom index is invalid!");
  }

  const AtomIndex index = _adjacencies.inner().addVertex(elementType);
  addBond(index, adjacentTo, bondType);
  // addBond handles the stereocenter update on adjacentTo

  _propagateGraphChange();

  return index;
}

void Molecule::Impl::addBond(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::addBond: A supplied index is invalid!");
  }

  if(a == b) {
    throw std::logic_error("Molecule::addBond: Cannot add a bond between identical indices!");
  }

  InnerGraph& inner = _adjacencies.inner();

  inner.addEdge(a, b, bondType);

  auto notifySubstituentAddition = [this](
    const AtomIndex toIndex,
    const AtomIndex addedIndex
  ) {
    /* Remove any BondStereocenters on adjacent edges of toIndex (no state
     * propagation possible yet TODO)
     */
    for(
      const BondIndex& adjacentEdge :
      boost::make_iterator_range(
        _adjacencies.bonds(toIndex)
      )
    ) {
      _stereocenters.try_remove(adjacentEdge);
    }

    if(auto atomStereocenterOption = _stereocenters.option(toIndex)) {
      // Re-rank around toIndex
      auto localRanking = rankPriority(toIndex);

      atomStereocenterOption->addSubstituent(
        _adjacencies,
        addedIndex,
        localRanking,
        determineLocalGeometry(toIndex, localRanking),
        Options::chiralStatePreservation
      );
    }
  };

  notifySubstituentAddition(a, b);
  notifySubstituentAddition(b, a);

  _propagateGraphChange();
}

void Molecule::Impl::assignStereocenter(
  const AtomIndex a,
  const boost::optional<unsigned>& assignmentOption
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::assignStereocenter: Supplied index is invalid!");
  }

  auto stereocenterOption = _stereocenters.option(a);

  if(!stereocenterOption) {
    throw std::out_of_range("assignStereocenter: No stereocenter at this index!");
  }

  if(
    assignmentOption
    && assignmentOption.value() >= stereocenterOption->numAssignments()
  ) {
    throw std::out_of_range("assignStereocenter: Invalid assignment index!");
  }

  stereocenterOption -> assign(assignmentOption);

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::assignStereocenter(
  const BondIndex& edge,
  const boost::optional<unsigned>& assignmentOption
) {
  if(!_isValidIndex(edge.first) || !_isValidIndex(edge.second)) {
    throw std::out_of_range("Molecule::assignStereocenter: Supplied bond atom indices is invalid!");
  }

  auto stereocenterOption = _stereocenters.option(edge);

  if(!stereocenterOption) {
    throw std::out_of_range("assignStereocenter: No stereocenter at this bond!");
  }

  if(
    assignmentOption
    && assignmentOption.value() >= stereocenterOption->numAssignments()
  ) {
    throw std::out_of_range("assignStereocenter: Invalid assignment index!");
  }

  stereocenterOption -> assign(assignmentOption);

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::assignStereocenterRandomly(const AtomIndex a) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::assignStereocenterRandomly: Supplied index is invalid!");
  }

  auto stereocenterOption = _stereocenters.option(a);

  if(!stereocenterOption) {
    throw std::out_of_range("assignStereocenterRandomly: No stereocenter at this index!");
  }

  stereocenterOption->assignRandom();

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::assignStereocenterRandomly(const BondIndex& e) {
  auto stereocenterOption = _stereocenters.option(e);

  if(!stereocenterOption) {
    throw std::out_of_range("assignStereocenterRandomly: No stereocenter at this edge!");
  }

  stereocenterOption->assignRandom();

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::removeAtom(const AtomIndex a) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::removeAtom: Supplied index is invalid!");
  }

  if(!graph().canRemove(a)) {
    throw std::logic_error("Removing this atom disconnects the graph!");
  }

  InnerGraph& inner = _adjacencies.inner();

  std::vector<AtomIndex> previouslyAdjacentVertices;
  auto adjacencyIterators = inner.adjacents(a);
  std::copy(
    adjacencyIterators.first,
    adjacencyIterators.second,
    std::back_inserter(previouslyAdjacentVertices)
  );

  // Remove all edges to and from this vertex
  inner.clearVertex(a);

  // Any stereocenter on this index must be dropped
  _stereocenters.try_remove(a);

  // Any adjacent bond stereocenters have to be dropped
  for(
    const BondIndex& adjacentEdge :
    boost::make_iterator_range(_adjacencies.bonds(a))
  ) {
    _stereocenters.try_remove(adjacentEdge);
  }

  // Remove the vertex itself
  inner.removeVertex(a);

  /* Removing the vertex invalidates some vertex descriptors, which are used
   * liberally in the stereocenter classes' state. We have to correct of all of
   * those to ensure that _propagateGraphChange works properly.
   */
  _stereocenters.propagateVertexRemoval(a);

  /* call removeSubstituent on all adjacent stereocenters, with
   * removalPlaceholder as the 'which' parameter, which is what
   * propagateVertexRemoval replaces the removed index with in the
   * stereocenters' internal state
   */
  for(const auto& indexToUpdate : previouslyAdjacentVertices) {
    if(auto stereocenterOption = _stereocenters.option(indexToUpdate)) {
      auto localRanking = rankPriority(indexToUpdate);

      // If index becomes terminal, drop the stereocenter
      if(localRanking.ligands.size() <= 1) {
        _stereocenters.remove(indexToUpdate);
        continue;
      }

      stereocenterOption->removeSubstituent(
        _adjacencies,
        InnerGraph::removalPlaceholder,
        localRanking,
        determineLocalGeometry(indexToUpdate, localRanking),
        Options::chiralStatePreservation
      );
    }

    // TODO BondStereocenter update
    /*if(_stereocenters.involving(indexToUpdate)) {
      if(_stereocenters.at(indexToUpdate) -> type() == Stereocenters::Type::AtomStereocenter) {
      } else {
        std::dynamic_pointer_cast<Stereocenters::BondStereocenter>(
          _stereocenters.at(indexToUpdate)
        ) -> removeSubstituent(
          indexToUpdate,
          Stereocenters::Stereocenter::removalPlaceholder
        );
      }
    }*/
  }

  _propagateGraphChange();
}

void Molecule::Impl::removeBond(
  const AtomIndex a,
  const AtomIndex b
) {
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::removeBond: Supplied index is invalid!");
  }

  InnerGraph& inner = _adjacencies.inner();

  auto edgeOption = inner.edgeOption(a, b);

  if(!edgeOption) {
    throw std::out_of_range("That bond does not exist!");
  }

  InnerGraph::Edge edgeToRemove = edgeOption.value();

  if(!inner.canRemove(edgeToRemove)) {
    throw std::logic_error("Removing this bond separates the molecule into two pieces!");
  }


  /* If there is an BondStereocenter on this edge, we have to drop it explicitly,
   * since _propagateGraphChange cannot iterate over a now-removed edge.
   */
  _stereocenters.try_remove(BondIndex {a, b});

  // Remove the edge
  inner.removeEdge(edgeToRemove);

  // Notify all immediately adjacent stereocenters of the removal
  auto notifyRemoval = [&](
    const AtomIndex indexToUpdate,
    const AtomIndex removedIndex
  ) {
    if(auto stereocenterOption = _stereocenters.option(indexToUpdate)) {
      auto localRanking = rankPriority(indexToUpdate);

      // In case the central atom becomes terminal, just drop the stereocenter
      if(localRanking.ligands.size() <= 1) {
        _stereocenters.remove(indexToUpdate);
        return;
      }

      stereocenterOption->removeSubstituent(
        _adjacencies,
        removedIndex,
        localRanking,
        determineLocalGeometry(indexToUpdate, localRanking),
        Options::chiralStatePreservation
      );
    }

    // TODO
    /*if(_stereocenters.involving(indexToUpdate)) {
      if(_stereocenters.at(indexToUpdate) -> type() == Stereocenters::Type::AtomStereocenter) {
      } else {
        std::dynamic_pointer_cast<Stereocenters::BondStereocenter>(
          _stereocenters.at(indexToUpdate)
        ) -> removeSubstituent(
          indexToUpdate,
          removedIndex
        );
      }
    }*/
  };

  notifyRemoval(a, b);
  notifyRemoval(b, a);

  /* All other cases, where there may be BondStereocenters or AtomStereocenters
   * on a or b, should be handled correctly by _propagateGraphChange.
   */

  _propagateGraphChange();
}

bool Molecule::Impl::setBondType(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::setBondType: A supplied index is invalid!");
  }

  if(bondType == BondType::Eta) {
    throw std::logic_error(
      "Do not manually change eta bond types, this dynamism is handled internally"
    );
  }

  InnerGraph& inner = _adjacencies.inner();

  auto edgeOption = inner.edgeOption(a, b);
  if(!edgeOption) {
    addBond(a, b, bondType);
    return false;
  }

  inner.bondType(edgeOption.value()) = bondType;
  _propagateGraphChange();
  return true;
}

void Molecule::Impl::setElementType(
  const AtomIndex a,
  const Delib::ElementType elementType
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::setElementType: This index is invalid!");
  }

  _adjacencies.inner().elementType(a) = elementType;
  _propagateGraphChange();
}

void Molecule::Impl::setGeometryAtAtom(
  const AtomIndex a,
  const Symmetry::Name symmetryName
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::setGeometryAtAtom: Supplied atom index is invalid");
  }

  auto stereocenterOption = _stereocenters.option(a);

  // If there is no stereocenter at this position yet, we have to create it
  if(!stereocenterOption) {
    RankingInformation localRanking = rankPriority(a);

    if(localRanking.ligands.size() != Symmetry::size(symmetryName)) {
      throw std::logic_error(
        "Molecule::setGeometryAtAtom: The size of the supplied symmetry is not "
        " the same as the number of determined ligands"
      );
    }

    // Add the stereocenter irrespective of how many assignments it has
    auto newStereocenter = AtomStereocenter {
      _adjacencies,
      symmetryName,
      a,
      localRanking
    };

    // Default-assign stereocenters with only one assignment
    if(newStereocenter.numAssignments() == 1) {
      newStereocenter.assign(0u);
    }

    _stereocenters.add(
      a,
      std::move(newStereocenter)
    );

    _propagateGraphChange();
    return;
  }

  if(
    Symmetry::size(stereocenterOption->getSymmetry())
    != Symmetry::size(symmetryName)
  ) {
    throw std::logic_error(
      "Molecule::setGeometryAtAtom: The size of the supplied symmetry is "
      "not the same as that of the existing stereocenter's current symmetry!"
    );
  }

  stereocenterOption->setSymmetry(symmetryName, _adjacencies);
  if(stereocenterOption->numStereopermutations() == 1) {
    assignStereocenter(a, 0);
    return;
  }
  _propagateGraphChange();
}

/* Information */
Symmetry::Name Molecule::Impl::determineLocalGeometry(
  const AtomIndex index,
  const RankingInformation& ranking
) const {
  if(!_isValidIndex(index)) {
    throw std::out_of_range("Molecule::determineLocalGeometry: Supplied index is invalid!");
  }

  if(graph().degree(index) <= 1) {
    throw std::logic_error(
      "Molecule::determineLocalGeometry: No geometries exist for terminal atoms"
    );
  }

  return LocalGeometry::determineLocalGeometry(
    _adjacencies,
    index,
    ranking
  );
}

std::string Molecule::Impl::dumpGraphviz() const {
  MolGraphWriter propertyWriter(&_adjacencies.inner());

  std::stringstream graphvizStream;

  boost::write_graphviz(
    graphvizStream,
    _adjacencies.inner().bgl(),
    propertyWriter,
    propertyWriter,
    propertyWriter
  );

  return graphvizStream.str();
}

const OuterGraph& Molecule::Impl::graph() const {
  return _adjacencies;
}

const StereocenterList& Molecule::Impl::stereocenters() const {
  return _stereocenters;
}

StereocenterList Molecule::Impl::inferStereocentersFromPositions(
  const AngstromWrapper& angstromWrapper,
  const boost::optional<
    std::vector<BondIndex>
  >& explicitBondStereocenterCandidatesOption
) const {
  const AtomIndex size = graph().N();
  StereocenterList stereocenters;

  Cycles cycleData = graph().cycles();

  for(unsigned candidateIndex = 0; candidateIndex < size; candidateIndex++) {
    RankingInformation localRanking = rankPriority(candidateIndex, {}, angstromWrapper);

    // Skip terminal atoms
    if(localRanking.ligands.size() <= 1) {
      continue;
    }

    const auto expectedGeometry = determineLocalGeometry(candidateIndex, localRanking);

    // Construct it
    auto stereocenter = AtomStereocenter {
      _adjacencies,
      expectedGeometry,
      candidateIndex,
      localRanking
    };

    pickyFit(
      stereocenter,
      _adjacencies,
      angstromWrapper,
      expectedGeometry
    );

    if(
      !disregardStereocenter(
        stereocenter,
        _adjacencies.elementType(candidateIndex),
        cycleData,
        Options::temperatureRegime
      )
    ) {
      stereocenters.add(
        candidateIndex,
        std::move(stereocenter)
      );
    }
  }

  const InnerGraph& inner = _adjacencies.inner();

  auto tryInstantiateBondStereocenter = [&](const InnerGraph::Edge& edgeIndex) -> void {
    InnerGraph::Vertex source = inner.source(edgeIndex),
                       target = inner.target(edgeIndex);

    auto sourceAtomStereocenterOption = stereocenters.option(source);
    auto targetAtomStereocenterOption = stereocenters.option(target);

    // There need to be assigned stereocenters on both vertices
    if(
      !sourceAtomStereocenterOption
      || !targetAtomStereocenterOption
      || sourceAtomStereocenterOption->assigned() == boost::none
      || targetAtomStereocenterOption->assigned() == boost::none
    ) {
      return;
    }

    const BondIndex bondIndex {source, target};

    // Construct a Stereocenter here
    auto newStereocenter = BondStereocenter {
      *sourceAtomStereocenterOption,
      *targetAtomStereocenterOption,
      bondIndex
    };

    newStereocenter.fit(
      angstromWrapper,
      *sourceAtomStereocenterOption,
      *targetAtomStereocenterOption
    );

    if(newStereocenter.assigned() != boost::none) {
      stereocenters.add(
        bondIndex,
        std::move(newStereocenter)
      );
    }
  };

  // Is there an explicit list of bonds on which to try BondStereocenter instantiation?
  if(explicitBondStereocenterCandidatesOption) {
    for(const BondIndex& bondIndex : *explicitBondStereocenterCandidatesOption) {
      // Test if the supplied edge exists first
      auto edge = inner.edgeOption(bondIndex.first, bondIndex.second);
      if(!edge) {
        throw std::out_of_range("Explicit bond stereocenter candidate edge does not exist!");
      }

      tryInstantiateBondStereocenter(*edge);
    }
  } else {
    // Every bond is a candidate
    for(
      const InnerGraph::Edge& edgeIndex :
      boost::make_iterator_range(inner.edges())
    ) {
      tryInstantiateBondStereocenter(edgeIndex);
    }
  }

  return stereocenters;
}

bool Molecule::Impl::modularCompare(
  const Molecule::Impl& other,
  const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
) const {
  const unsigned thisNumAtoms = graph().N();

  if(thisNumAtoms != other.graph().N()) {
    return false;
  }

  // Compare number of bonds too
  if(graph().B() != other.graph().B()) {
    return false;
  }

  /* boost isomorphism will allocate a vector of size maxHash, this is dangerous
   * as the maximum hash can be immense, another post-processing step is needed
   * for the calculated hashes to decrease the spatial requirements
   *
   * This maps the hashes to an incremented number
   */
  std::vector<hashes::HashType> thisHashes, otherHashes;
  hashes::HashType maxHash;

  std::tie(thisHashes, otherHashes, maxHash) = hashes::narrow(
    hashes::generate(graph().inner(), stereocenters(), comparisonBitmask),
    hashes::generate(other.graph().inner(), other.stereocenters(), comparisonBitmask)
  );

  // Where the corresponding index from the other graph is stored
  std::vector<AtomIndex> indexMap(thisNumAtoms);

  const auto& thisBGL = _adjacencies.inner().bgl();
  const auto& otherBGL = other._adjacencies.inner().bgl();

  bool isomorphic = boost::isomorphism(
    thisBGL,
    otherBGL,
    boost::make_safe_iterator_property_map(
      indexMap.begin(),
      thisNumAtoms,
      boost::get(boost::vertex_index, thisBGL)
    ),
    hashes::LookupFunctor(thisHashes),
    hashes::LookupFunctor(otherHashes),
    maxHash,
    boost::get(boost::vertex_index, thisBGL),
    boost::get(boost::vertex_index, otherBGL)
  );

  return isomorphic;
}

RankingInformation Molecule::Impl::rankPriority(
  const AtomIndex a,
  const std::vector<AtomIndex>& excludeAdjacent,
  const boost::optional<AngstromWrapper>& positionsOption
) const {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Supplied atom index is invalid!");
  }

  RankingInformation rankingResult;

  // Expects that bond types are set properly, complains otherwise
  rankingResult.ligands = GraphAlgorithms::ligandSiteGroups(
    _adjacencies.inner(),
    a,
    excludeAdjacent
  );

  std::string molGraphviz;
#ifndef NDEBUG
  molGraphviz = dumpGraphviz();
#endif

  // Rank the substituents
  auto expandedTree = RankingTree(
    graph(),
    graph().cycles(),
    stereocenters(),
    molGraphviz,
    a,
    excludeAdjacent,
    RankingTree::ExpansionOption::OnlyRequiredBranches,
    positionsOption
  );

  rankingResult.sortedSubstituents = expandedTree.getRanked();

  rankingResult.ligandsRanking = RankingInformation::rankLigands(
    rankingResult.ligands,
    rankingResult.sortedSubstituents
  );

  // Find links between them
  rankingResult.links = GraphAlgorithms::substituentLinks(
    _adjacencies.inner(),
    graph().cycles(),
    a,
    rankingResult.ligands,
    excludeAdjacent
  );

  return rankingResult;
}

bool Molecule::Impl::operator == (const Impl& other) const {
  // Operator == performs the most strict comparison possible
  return modularCompare(
    other,
    temple::make_bitmask(AtomEnvironmentComponents::ElementTypes)
      | AtomEnvironmentComponents::BondOrders
      | AtomEnvironmentComponents::Symmetries
      | AtomEnvironmentComponents::Stereopermutations
  );
}

bool Molecule::Impl::operator != (const Impl& other) const {
  return !(*this == other);
}

} // namespace molassembler
