#include "molassembler/Molecule.h"

#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/graph_utility.hpp"
#include "chemical_symmetries/ConstexprProperties.h"
#include "Delib/Constants.h"
#include "Delib/ElementTypeCollection.h"
#include "temple/constexpr/Numeric.h"
#include "temple/Containers.h"
#include "temple/Optionals.h"

#include "molassembler/Cycles.h"
#include "molassembler/Graph/InnerGraph.h"
#include "molassembler/Graph/GraphAlgorithms.h"
#include "molassembler/Modeling/CommonTrig.h"
#include "molassembler/Modeling/LocalGeometryModel.h"
#include "molassembler/Molecule/MolGraphWriter.h"
#include "molassembler/Molecule/RankingTree.h"
#include "molassembler/RankingInformation.h"
#include "molassembler/StereocenterList.h"
#include "molassembler/Options.h"
#include "molassembler/OuterGraph.h"

namespace molassembler {

/* Impl definition */
struct Molecule::Impl {
  OuterGraph _adjacencies;
  StereocenterList _stereocenters;

/* "Private" helpers */
  //! Generates a list of stereocenters based on graph properties alone
  StereocenterList _detectStereocenters() const;

  //! Ensures basic expectations about what constitutes a Molecule are met
  void _ensureModelInvariants() const;

  //! Returns whether the specified index is valid or not
  bool _isValidIndex(AtomIndex index) const;

  //! Returns whether an edge is double, aromtic, triple or higher bond order
  bool _isGraphBasedBondStereocenterCandidate(const InnerGraph::Edge& e) const;

  //! Updates the molecule's StereocenterList after a graph modification
  void _propagateGraphChange();


//!@name Constructors
//!@{
  //! Default constructor
  Impl() noexcept;

  //! Diatomic constructor
  Impl(
    Delib::ElementType a,
    Delib::ElementType b,
    BondType bondType
  ) noexcept;

  //! Graph-only constructor
  explicit Impl(OuterGraph graph);

  //! Graph and positions constructor
  Impl(
    OuterGraph graph,
    const AngstromWrapper& positions
  );

  //! Graph and stereocenters constructor
  Impl(
    OuterGraph graph,
    StereocenterList stereocenters
  );
//!@}

//!@name Modifiers
//!@{
  //! Adds an atom by attaching it to an existing atom.
  AtomIndex addAtom(
    Delib::ElementType elementType,
    AtomIndex adjacentTo,
    BondType bondType
  );

  //! Adds a bond between existing atoms.
  void addBond(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  /*! Sets the stereocenter assignment at a particular atom
   *
   * This sets the stereocenter assignment at a specific atom index. For this,
   * a stereocenter must be instantiated and contained in the StereocenterList
   * returned by stereocenters(). The supplied assignment must be either
   * boost::none or smaller than stereocenterPtr->numAssignments().
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenter(
    AtomIndex a,
    const boost::optional<unsigned>& assignment
  );

  /*! Assigns a stereocenter stereopermutation at random
   *
   * This sets the stereocetner assignment at a specific index, taking relative
   * statistical occurence weights of each stereopermutation into account.
   *
   * \pre There must be an AtomStereocenter at the passed index
   *
   * \throws If no AtomStereocenter exists at the passed index
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators to its members invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenterRandomly(AtomIndex a);

  /*! Assigns a bond stereocenter to a random assignment
   *
   * \pre There must be a BondStereocenter at the passed edge
   * \throws If no BondStereocenter exists at the passed edge
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators to its members invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenterRandomly(const BondIndex& e);

  /*! Removes an atom from the graph, including bonds to it.
   *
   * Removes an atom from the molecular graph, including bonds to the atom,
   * after checking that removing it is safe, i.e. the removal does not
   * disconnect the graph.
   *
   * \throws if the supplied index is invalid or isSafeToRemoveAtom returns false.
   */
  void removeAtom(AtomIndex a);

  /*!
   * Removes an atom after checking if removing that bond is safe, i.e. does not
   * disconnect the graph. An example of bonds that can always be removed are
   * ring-closing bonds, since they never disconnect the molecular graph.
   *
   * \throws if isSafeToRemoveBond returns false.
   *
   * \note It is not safe to remove a bond just because one of the involved
   * atoms is terminal, since that atom would then be disconnected from the
   * rest of the molecule. This function merely removes a bond from the graph.
   * It is, however, considered safe to remove the terminal vertex, which
   * involves removing the bond to it.
   */
  void removeBond(AtomIndex a, AtomIndex b);

  //! Changes an existing bond's type
  bool setBondType(
    AtomIndex a,
    AtomIndex b,
    BondType bondType
  );

  //! Changes an existing atom's element type
  void setElementType(
    AtomIndex a,
    Delib::ElementType elementType
  );

  /*! Sets the local geometry at an atom index
   *
   * This sets the local geometry at a specific atom index. There are a number
   * of cases that this function treats differently, besides faulty arguments:
   * If there is already a AtomStereocenter instantiated at this atom index, its
   * underlying symmetry is altered. If there is no AtomStereocenter at
   * this index, one is instantiated. In all cases, new or modified
   * stereocenters are default-assigned if there is only one possible
   * assignment.
   * \throws if
   *   - the supplied atomic index is invalid
   *   - there is an BondStereocenter at that index
   *   - or the provided symmetry is a different size than that of an existing
   *     AtomStereocenter or the expected symmetry
   */
  void setGeometryAtAtom(
    AtomIndex a,
    Symmetry::Name symmetryName
  );
//!@}

//!@name Information
//!@{
  /*! Determines what the local geometry at a non-terminal atom ought to be
   *
   * Returns the expected symmetry name at a non-terminal atom.
   * \throws if the supplied atomic index is invalid
   */
  Symmetry::Name determineLocalGeometry(
    AtomIndex index,
    const RankingInformation& ranking
  ) const;

  //! Returns a graphivz string representation of the molecule
  std::string dumpGraphviz() const;

  Cycles cycles() const;

  //! Provides read-only access to the graph member
  const OuterGraph& graph() const;

  //! Provides read-only access to the list of stereocenters
  const StereocenterList& stereocenters() const;

  //! Returns the number of adjacencies of an atomic position
  unsigned getNumAdjacencies(AtomIndex a) const;

  StereocenterList inferStereocentersFromPositions(const AngstromWrapper& angstromWrapper) const;

  //! Modular comparison of this Impl with another.
  bool modularCompare(
    const Impl& other,
    const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
  ) const;

  RankingInformation rankPriority(
    AtomIndex a,
    const std::vector<AtomIndex>& excludeAdjacent = {},
    const boost::optional<AngstromWrapper>& positionsOption = boost::none
  ) const;
//!@}

//!@name Operators
//!@{
  //! Equality operator, performs most strict equality comparison
  bool operator == (const Impl& other) const;
  bool operator != (const Impl& other) const;
//!@}
};

/* Impl members */
/* Private members */
StereocenterList Molecule::Impl::_detectStereocenters() const {
  StereocenterList stereocenterList;

  Cycles cycleData = cycles();
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
    || bondType == BondType::Aromatic
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

  Cycles cycleData = cycles();

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

      // Unconditionally propagate the graph change
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
  const AngstromWrapper& positions
) : _adjacencies(std::move(graph))
{
  GraphAlgorithms::findAndSetEtaBonds(_adjacencies.inner());
  _stereocenters = inferStereocentersFromPositions(positions);
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
  const boost::optional<unsigned>& assignment
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::assignStereocenter: Supplied index is invalid!");
  }

  auto stereocenterOption = _stereocenters.option(a);

  if(!stereocenterOption) {
    throw std::logic_error("assignStereocenter: No stereocenter at this index!");
  }

  if(assignment >= stereocenterOption->numAssignments()) {
    throw std::logic_error("assignStereocenter: Invalid assignment index!");
  }

  stereocenterOption -> assign(assignment);

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::assignStereocenterRandomly(const AtomIndex a) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::assignStereocenterRandomly: Supplied index is invalid!");
  }

  auto stereocenterOption = _stereocenters.option(a);

  if(!stereocenterOption) {
    throw std::logic_error("assignStereocenterRandomly: No stereocenter at this index!");
  }

  stereocenterOption->assignRandom();

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::assignStereocenterRandomly(const BondIndex& e) {
  auto stereocenterOption = _stereocenters.option(e);

  if(!stereocenterOption) {
    throw std::logic_error("assignStereocenterRandomly: No stereocenter at this edge!");
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
    throw std::logic_error("That bond does not exist!");
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

      // In case the CNS central atom becomes terminal, just drop the stereocenter
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
    const Symmetry::Name expectedSymmetry = determineLocalGeometry(a, localRanking);

    if(localRanking.ligands.size() != Symmetry::size(symmetryName)) {
      throw std::logic_error(
        "Molecule::setGeometryAtAtom: The size of the supplied symmetry is not "
        " the same as the number of determined ligands"
      );
    }

    /* In case the expected symmetry is the same as the supplied symmetry, this
     * is a no-op, since it adds no information to the Molecule (this
     * stereocenter would otherwise be created during DG initialization)
     */
    if(expectedSymmetry != symmetryName) {
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
    }
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

  if(getNumAdjacencies(index) <= 1) {
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

Cycles Molecule::Impl::cycles() const {
  return Cycles {_adjacencies};
}

const OuterGraph& Molecule::Impl::graph() const {
  return _adjacencies;
}

const StereocenterList& Molecule::Impl::stereocenters() const {
  return _stereocenters;
}

unsigned Molecule::Impl::getNumAdjacencies(const AtomIndex a) const {
  return _adjacencies.inner().degree(a);
}

StereocenterList Molecule::Impl::inferStereocentersFromPositions(
  const AngstromWrapper& angstromWrapper
) const {
  const AtomIndex size = graph().N();
  StereocenterList stereocenters;

  Cycles cycleData = cycles();

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

  for(
    const InnerGraph::Edge& edgeIndex :
    boost::make_iterator_range(inner.edges())
  ) {
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
      continue;
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

  auto thisHashes = hashes::generate(graph().inner(), stereocenters(), comparisonBitmask);
  auto otherHashes = hashes::generate(other.graph().inner(), other.stereocenters(), comparisonBitmask);

  /* boost isomorphism will allocate a vector of size maxHash, this is dangerous
   * as the maximum hash can be immense, another post-processing step is needed
   * for the calculated hashes to decrease the spatial requirements
   *
   * This maps the hashes to an incremented number
   */
  auto maxHash = hashes::regularize(thisHashes, otherHashes);

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
    cycles(),
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
    cycles(),
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


/* Molecule interface to Impl call forwards */
Molecule::Molecule() noexcept : _pImpl(
  std::make_unique<Impl>()
) {}

Molecule::Molecule(Molecule&& other) noexcept = default;
Molecule& Molecule::operator = (Molecule&& rhs) noexcept = default;

Molecule::Molecule(const Molecule& other) : _pImpl(
  std::make_unique<Impl>(*other._pImpl)
) {}
Molecule& Molecule::operator = (const Molecule& rhs) {
  *_pImpl = *rhs._pImpl;
  return *this;
}

Molecule::~Molecule() = default;

Molecule::Molecule(
  const Delib::ElementType a,
  const Delib::ElementType b,
  const BondType bondType
) noexcept : _pImpl(
  std::make_unique<Impl>(a, b, bondType)
) {}

Molecule::Molecule(OuterGraph graph) : _pImpl(
  std::make_unique<Impl>(std::move(graph))
) {}

Molecule::Molecule(
  OuterGraph graph,
  const AngstromWrapper& positions
) : _pImpl(
  std::make_unique<Impl>(
    std::move(graph),
    positions
  )
) {}

Molecule::Molecule(
  OuterGraph graph,
  StereocenterList stereocenters
) : _pImpl(
  std::make_unique<Impl>(
    std::move(graph),
    std::move(stereocenters)
  )
) {}

/* Modifiers */
AtomIndex Molecule::addAtom(
  const Delib::ElementType elementType,
  const AtomIndex adjacentTo,
  const BondType bondType
) {
  return _pImpl->addAtom(elementType, adjacentTo, bondType);
}

void Molecule::addBond(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  _pImpl->addBond(a, b, bondType);
}

void Molecule::assignStereocenter(
  const AtomIndex a,
  const boost::optional<unsigned>& assignment
) {
  _pImpl->assignStereocenter(a, assignment);
}

void Molecule::assignStereocenterRandomly(const AtomIndex a) {
  _pImpl->assignStereocenterRandomly(a);
}

void Molecule::assignStereocenterRandomly(const BondIndex& e) {
  _pImpl->assignStereocenterRandomly(e);
}

void Molecule::removeAtom(const AtomIndex a) {
  _pImpl->removeAtom(a);
}

void Molecule::removeBond(
  const AtomIndex a,
  const AtomIndex b
) {
  _pImpl->removeBond(a, b);
}

bool Molecule::setBondType(
  const AtomIndex a,
  const AtomIndex b,
  const BondType bondType
) {
  return _pImpl->setBondType(a, b, bondType);
}

void Molecule::setElementType(
  const AtomIndex a,
  const Delib::ElementType elementType
) {
  _pImpl->setElementType(a, elementType);
}

void Molecule::setGeometryAtAtom(
  const AtomIndex a,
  const Symmetry::Name symmetryName
) {
  _pImpl->setGeometryAtAtom(a, symmetryName);
}


/* Information */
Symmetry::Name Molecule::determineLocalGeometry(
  const AtomIndex index,
  const RankingInformation& ranking
) const {
  return _pImpl->determineLocalGeometry(index, ranking);
}

std::string Molecule::dumpGraphviz() const {
  return _pImpl->dumpGraphviz();
}

Cycles Molecule::cycles() const {
  return _pImpl->cycles();
}

const OuterGraph& Molecule::graph() const {
  return _pImpl->graph();
}

const StereocenterList& Molecule::stereocenters() const {
  return _pImpl->stereocenters();
}

unsigned Molecule::getNumAdjacencies(const AtomIndex a) const {
  return _pImpl->getNumAdjacencies(a);
}

StereocenterList Molecule::inferStereocentersFromPositions(
  const AngstromWrapper& angstromWrapper
) const {
  return _pImpl->inferStereocentersFromPositions(angstromWrapper);
}

bool Molecule::modularCompare(
  const Molecule& other,
  const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
) const {
  return _pImpl->modularCompare(*other._pImpl, comparisonBitmask);
}

RankingInformation Molecule::rankPriority(
  const AtomIndex a,
  const std::vector<AtomIndex>& excludeAdjacent,
  const boost::optional<AngstromWrapper>& positionsOption
) const {
  return _pImpl->rankPriority(a, excludeAdjacent, positionsOption);
}

/* Operators */
bool Molecule::operator == (const Molecule& other) const {
  return *_pImpl == *other._pImpl;
}

bool Molecule::operator != (const Molecule& other) const {
  return *_pImpl != *other._pImpl;
}


} // namespace molassembler

std::ostream& operator << (
  std::ostream& os,
  const molassembler::Molecule& molecule
) {
  const auto& stereocenters = molecule.stereocenters();

  if(!stereocenters.empty()) {
    os << "Stereocenter information:\n";

    for(const auto& stereocenter : stereocenters.atomStereocenters()) {
      os << stereocenter.info() << "\n";
    }

    for(const auto& stereocenter : stereocenters.bondStereocenters()) {
      os << stereocenter.info() << "\n";
    }
  }

  return os;
}
