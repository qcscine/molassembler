#include "temple/Optionals.h"

#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/graph_utility.hpp"

#include "Delib/Constants.h"
#include "Delib/ElementTypeCollection.h"

#include "chemical_symmetries/ConstexprProperties.h"

#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"

#include "detail/CommonTrig.h"
#include "detail/MolGraphWriter.h"
#include "Cycles.h"
#include "GraphAlgorithms.h"
#include "GraphHelpers.h"
#include "LocalGeometryModel.h"
#include "Log.h"
#include "Molecule.h"
#include "RankingInformation.h"
#include "RankingTree.h"
#include "StereocenterList.h"
#include "Options.h"

namespace molassembler {

/* Impl definition */
struct Molecule::Impl {
  GraphType _adjacencies;
  StereocenterList _stereocenters;

/* "Private" helpers */
  //! Generates a list of stereocenters based on graph properties alone
  StereocenterList _detectStereocenters() const;

  //! Ensures basic expectations about what constitutes a Molecule are met
  void _ensureModelInvariants() const;

  //! Returns whether the specified index is valid or not
  bool _isValidIndex(AtomIndexType index) const;

  //! Returns whether an edge is double, aromtic, triple or higher bond order
  bool _isGraphBasedBondStereocenterCandidate(const GraphType::edge_descriptor& e) const;

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
  explicit Impl(GraphType graph);

  //! Graph and positions constructor
  Impl(
    GraphType graph,
    const AngstromWrapper& positions
  );

  //! Graph and stereocenters constructor
  Impl(
    GraphType graph,
    StereocenterList stereocenters
  );
//!@}

//!@name Modifiers
//!@{
  //! Adds an atom by attaching it to an existing atom.
  AtomIndexType addAtom(
    Delib::ElementType elementType,
    AtomIndexType adjacentTo,
    BondType bondType
  );

  //! Adds a bond between existing atoms.
  void addBond(
    AtomIndexType a,
    AtomIndexType b,
    BondType bondType
  );

  /*! Sets the stereocenter assignment at a particular atom
   *
   * This sets the stereocenter assignment at a specific atom index. For this,
   * a stereocenter must be instantiated and contained in the StereocenterList
   * returned by getStereocenterList(). The supplied assignment must be either
   * boost::none or smaller than stereocenterPtr->numAssignments().
   *
   * \note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenter(
    AtomIndexType a,
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
  void assignStereocenterRandomly(AtomIndexType a);

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
  void assignStereocenterRandomly(const GraphType::edge_descriptor& e);

  /*! Removes an atom from the graph, including bonds to it.
   *
   * Removes an atom from the molecular graph, including bonds to the atom,
   * after checking that removing it is safe, i.e. the removal does not
   * disconnect the graph.
   *
   * \throws if the supplied index is invalid or isSafeToRemoveAtom returns false.
   */
  void removeAtom(AtomIndexType a);

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
  void removeBond(AtomIndexType a, AtomIndexType b);

  //! Changes an existing bond's type
  bool setBondType(
    AtomIndexType a,
    AtomIndexType b,
    BondType bondType
  );

  //! Changes an existing atom's element type
  void setElementType(
    AtomIndexType a,
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
    AtomIndexType a,
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
    AtomIndexType index,
    const RankingInformation& ranking
  ) const;

  //! Returns a graphivz string representation of the molecule
  std::string dumpGraphviz() const;

  //! Get an edge descriptor for a pair of indices
  GraphType::edge_descriptor edge(
    AtomIndexType a,
    AtomIndexType b
  ) const;

  /*! Fetches the atomic indices of vertices adjacent to a particular index
   *
   * Fetches the atomic indices of vertices adjacent to a particular index.
   * \throws if the supplied atomic index is invalid
   */
  std::vector<AtomIndexType> getAdjacencies(AtomIndexType a) const;

  //! Fetches the optional bond type between two atom indices
  boost::optional<BondType> getBondType(
    AtomIndexType a,
    AtomIndexType b
  ) const;

  BondType getBondType(const GraphType::edge_descriptor& edge) const;

  Cycles getCycleData() const;

  /*! Returns the element type of an atomic index
   *
   * Returns the element type of an atomic index
   * \throws if the atomic index is invalid
   */
  Delib::ElementType getElementType(AtomIndexType index) const;

  //! Returns a collection detailing all element types
  Delib::ElementTypeCollection getElementCollection() const;

  //! Provides read-only access to the graph member
  const GraphType& getGraph() const;

  //! Provides read-only access to the list of stereocenters
  const StereocenterList& getStereocenterList() const;

  //! Returns the number of adjacencies of an atomic position
  unsigned getNumAdjacencies(AtomIndexType a) const;

  StereocenterList inferStereocentersFromPositions(const AngstromWrapper& angstromWrapper) const;

  //! Checks if two atomic indices are connected by a bond
  bool isAdjacent(
    AtomIndexType a,
    AtomIndexType b
  ) const;

  //! An atom is considered removable if it isn't an articulation vertex
  bool isSafeToRemoveAtom(AtomIndexType a) const;

  //! A bond is considered removable if it isn't a bridge edge
  bool isSafeToRemoveBond(AtomIndexType a, AtomIndexType b) const;

  //! Modular comparison of this Impl with another.
  bool modularCompare(
    const Impl& other,
    const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
  ) const;

  //! Fetch the number of atoms
  unsigned numAtoms() const;

  //! Fetch the number of bonds
  unsigned numBonds() const;

  RankingInformation rankPriority(
    AtomIndexType a,
    const std::set<AtomIndexType>& excludeAdjacent = {},
    const boost::optional<AngstromWrapper>& positionsOption = boost::none
  ) const;

  //! Get the vertex indices on both ends of a graph edge
  std::array<AtomIndexType, 2> vertices(const GraphType::edge_descriptor& edge) const;
//!@}

//!@name Iterators
//!@{
  //! Returns a range-for temporary object iterating through all atom indices
  RangeForTemporary<GraphType::vertex_iterator> iterateAtoms() const;

  /*! Returns a range-for temporary object allowing c++11 style for loop
   * iteration through an atom's adjacencies
   */
  RangeForTemporary<GraphType::adjacency_iterator> iterateAdjacencies(
    AtomIndexType a
  ) const;

  /*! Returns a range-for temporary object allowing c++11-style for loop
   * iteration through edges
   */
  RangeForTemporary<GraphType::edge_iterator> iterateEdges() const;

  /*! Returns a range-for temporary object allowing c++11-style for loop
   * iteration through edges around a specific atom
   */
  RangeForTemporary<GraphType::out_edge_iterator> iterateEdges(AtomIndexType a) const;
//!@}


//!@name Operators
//!@{
  //! Returns the adjacencies of the specified atom index
  RangeForTemporary<GraphType::adjacency_iterator> operator [] (
    AtomIndexType a
  ) const;

  //! Equality operator, performs most strict equality comparison
  bool operator == (const Impl& other) const;
  bool operator != (const Impl& other) const;
//!@}

};

/* Impl members */
/* Private members */
StereocenterList Molecule::Impl::_detectStereocenters() const {
  StereocenterList stereocenterList;

  Cycles cycleData = getCycleData();
  // Find AtomStereocenters
  for(
    AtomIndexType candidateIndex = 0;
    candidateIndex < numAtoms();
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
        getElementType(candidateIndex),
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

  // Find BondStereocenters
  for(const auto& edgeIndex : graph::edges(_adjacencies)) {
    // Skip edges that cannot be candidates
    if(!_isGraphBasedBondStereocenterCandidate(edgeIndex)) {
      continue;
    }

    AtomIndexType source = boost::source(edgeIndex, _adjacencies),
                  target = boost::target(edgeIndex, _adjacencies);

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

    // Construct a Stereocenter here
    auto newStereocenter = BondStereocenter {
      *sourceAtomStereocenterOption,
      *targetAtomStereocenterOption,
      edgeIndex
    };

    if(newStereocenter.numAssignments() > 1) {
      stereocenterList.add(
        edgeIndex,
        std::move(newStereocenter)
      );
    }
  }

  return stereocenterList;
}

void Molecule::Impl::_ensureModelInvariants() const {
  if(GraphAlgorithms::numConnectedComponents(_adjacencies) > 1) {
    throw std::logic_error("Molecules must be a single connected component. The supplied graph has multiple");
  }

  if(numAtoms() < 2) {
    throw std::logic_error("Molecules must consist of at least two atoms!");
  }
}

bool Molecule::Impl::_isValidIndex(const AtomIndexType index) const {
  return index < numAtoms();
}

bool Molecule::Impl::_isGraphBasedBondStereocenterCandidate(
  const GraphType::edge_descriptor& e
) const {
  return (
    _adjacencies[e].bondType == BondType::Double
    || _adjacencies[e].bondType == BondType::Aromatic
    || _adjacencies[e].bondType == BondType::Triple
    || _adjacencies[e].bondType == BondType::Quadruple
    || _adjacencies[e].bondType == BondType::Quintuple
    || _adjacencies[e].bondType == BondType::Sextuple
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

  GraphAlgorithms::findAndSetEtaBonds(_adjacencies);

  /* TODO
   * - Need state propagation for BondStereocenters, anything else is madness
   */

  Cycles cycleData = getCycleData();

  for(const auto vertex : iterateAtoms()) {
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
          getElementType(vertex),
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
          getElementType(vertex),
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
  graph::addAtom(a, _adjacencies);
  graph::addAtom(b, _adjacencies);

  // Although addBond is potentially-throwing, it never will here
  addBond(0, 1, bondType);
}

Molecule::Impl::Impl(GraphType graph) {
  // boost::adjacency_list doesn't implement moves, so use the swap member
  _adjacencies.swap(graph);

  // Initialization
  GraphAlgorithms::findAndSetEtaBonds(_adjacencies);
  _stereocenters = _detectStereocenters();
  _ensureModelInvariants();
}

Molecule::Impl::Impl(
  GraphType graph,
  const AngstromWrapper& positions
) {
  // boost::adjacency_list doesn't implement moves, so use the swap member
  _adjacencies.swap(graph);

  GraphAlgorithms::findAndSetEtaBonds(_adjacencies);
  _stereocenters = inferStereocentersFromPositions(positions);
  _ensureModelInvariants();
}

Molecule::Impl::Impl(
  GraphType graph,
  StereocenterList stereocenters
) : _stereocenters(std::move(stereocenters)) {
  // boost::adjacency_list doesn't implement moves, so use the swap member
  _adjacencies.swap(graph);

  // Initialization
  _ensureModelInvariants();
}

/* Modifiers */
AtomIndexType Molecule::Impl::addAtom(
  const Delib::ElementType elementType,
  const AtomIndexType adjacentTo,
  const BondType bondType
) {
  if(!_isValidIndex(adjacentTo)) {
    throw std::out_of_range("Molecule::addAtom: Supplied atom index is invalid!");
  }

  const auto index = graph::addAtom(elementType, _adjacencies);
  addBond(index, adjacentTo, bondType);
  // addBond handles the stereocenter update on adjacentTo

  _propagateGraphChange();

  return index;
}

void Molecule::Impl::addBond(
  const AtomIndexType a,
  const AtomIndexType b,
  const BondType bondType
) {
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::addBond: A supplied index is invalid!");
  }

  if(a == b) {
    throw std::logic_error("Molecule::addBond: Cannot add a bond between identical indices!");
  }

  graph::addBond(a, b, bondType, _adjacencies);

  auto notifySubstituentAddition = [this](
    const AtomIndexType toIndex,
    const AtomIndexType addedIndex
  ) {
    /* Remove any BondStereocenters on adjacent edges of toIndex (no state
     * propagation possible yet TODO)
     */
    for(const auto& adjacentEdge : graph::edges(toIndex, _adjacencies)) {
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
  const AtomIndexType a,
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

void Molecule::Impl::assignStereocenterRandomly(const AtomIndexType a) {
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

void Molecule::Impl::assignStereocenterRandomly(const GraphType::edge_descriptor& e) {
  auto stereocenterOption = _stereocenters.option(e);

  if(!stereocenterOption) {
    throw std::logic_error("assignStereocenterRandomly: No stereocenter at this edge!");
  }

  stereocenterOption->assignRandom();

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::removeAtom(const AtomIndexType a) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::removeAtom: Supplied index is invalid!");
  }

  if(!isSafeToRemoveAtom(a)) {
    throw std::logic_error("Removing this atom disconnects the graph!");
  }

  auto previouslyAdjacentVertices = getAdjacencies(a);

  // Remove all edges to and from this vertex
  boost::clear_vertex(a, _adjacencies);

  // Any stereocenter on this index must be dropped
  _stereocenters.try_remove(a);

  // Any adjacent bond stereocenters have to be dropped
  for(const auto& adjacentEdge : graph::edges(a, _adjacencies)) {
    _stereocenters.try_remove(adjacentEdge);
  }

  // Remove the vertex itself
  boost::remove_vertex(a, _adjacencies);

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
        Stereocenters::removalPlaceholder,
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
  const AtomIndexType a,
  const AtomIndexType b
) {
  if(!_isValidIndex(a) || !_isValidIndex(b)) {
    throw std::out_of_range("Molecule::removeBond: Supplied index is invalid!");
  }

  if(!isSafeToRemoveBond(a, b)) {
    throw std::logic_error("Removing this bond separates the molecule into two pieces!");
  }

  // Find edge
  auto edgePair = boost::edge(a, b, _adjacencies);

  if(!edgePair.second) {
    throw std::logic_error("That bond does not exist!");
  }

  /* If there is an BondStereocenter on this edge, we have to drop it explicitly,
   * since _propagateGraphChange cannot iterate over a now-removed edge.
   */
  _stereocenters.try_remove(edgePair.first);

  // Remove the edge
  boost::remove_edge(edgePair.first, _adjacencies);

  // Notify all immediately adjacent stereocenters of the removal
  auto notifyRemoval = [&](
    const AtomIndexType indexToUpdate,
    const AtomIndexType removedIndex
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
  const AtomIndexType a,
  const AtomIndexType b,
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

  auto edgePair = boost::edge(a, b, _adjacencies);
  if(edgePair.second) {
    _adjacencies[edgePair.first].bondType = bondType;
  } else {
    addBond(a, b, bondType);
  }

  _propagateGraphChange();

  return edgePair.second;
}

void Molecule::Impl::setElementType(
  const AtomIndexType a,
  const Delib::ElementType elementType
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::setElementType: This index is invalid!");
  }

  _adjacencies[a].elementType = elementType;
  _propagateGraphChange();
}

void Molecule::Impl::setGeometryAtAtom(
  const AtomIndexType a,
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
  const AtomIndexType index,
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

//! Get an edge descriptor for a pair of indices
GraphType::edge_descriptor Molecule::Impl::edge(
  const AtomIndexType a,
  const AtomIndexType b
) const {
  return graph::edge(a, b, _adjacencies);
}

std::vector<AtomIndexType> Molecule::Impl::getAdjacencies(const AtomIndexType a) const {
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

boost::optional<BondType> Molecule::Impl::getBondType(
  const AtomIndexType a,
  const AtomIndexType b
) const {
  if(auto edgeOption = graph::edgeOption(a, b, _adjacencies)) {
    return graph::bondType(edgeOption.value(), _adjacencies);
  }

  return boost::none;
}

BondType Molecule::Impl::getBondType(const GraphType::edge_descriptor& edge) const {
  return graph::bondType(edge, _adjacencies);
}

Cycles Molecule::Impl::getCycleData() const {
  return Cycles {_adjacencies};
}

Delib::ElementType Molecule::Impl::getElementType(const AtomIndexType index) const {
  if(!_isValidIndex(index)) {
    throw std::out_of_range("Molecule::getElementType: Supplied index is invalid");
  }

  return graph::elementType(index, _adjacencies);
}

Delib::ElementTypeCollection Molecule::Impl::getElementCollection() const {
  AtomIndexType N = numAtoms();
  Delib::ElementTypeCollection elements;
  elements.reserve(N);

  for(AtomIndexType i = 0; i < N; ++i) {
    elements.push_back(
      graph::elementType(i, _adjacencies)
    );
  }

  return elements;
}

const GraphType& Molecule::Impl::getGraph() const {
  return _adjacencies;
}

const StereocenterList& Molecule::Impl::getStereocenterList() const {
  return _stereocenters;
}

unsigned Molecule::Impl::getNumAdjacencies(const AtomIndexType a) const {
  return graph::numAdjacencies(a, _adjacencies);
}

StereocenterList Molecule::Impl::inferStereocentersFromPositions(
  const AngstromWrapper& angstromWrapper
) const {
  StereocenterList stereocenters;

  Cycles cycleData = getCycleData();

  for(unsigned candidateIndex = 0; candidateIndex < numAtoms(); candidateIndex++) {
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
        graph::elementType(candidateIndex, _adjacencies),
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

  for(const auto& edgeIndex : graph::edges(_adjacencies)) {
    auto source = graph::source(edgeIndex, _adjacencies),
         target = graph::target(edgeIndex, _adjacencies);

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

    // Construct a Stereocenter here
    auto newStereocenter = BondStereocenter {
      *sourceAtomStereocenterOption,
      *targetAtomStereocenterOption,
      edgeIndex
    };

    std::cout << "Fitting on " << source << " - " << target << "\n";
    newStereocenter.fit(
      angstromWrapper,
      *sourceAtomStereocenterOption,
      *targetAtomStereocenterOption
    );

    if(newStereocenter.assigned() != boost::none) {
      stereocenters.add(
        edgeIndex,
        std::move(newStereocenter)
      );
    }
  }

  return stereocenters;
}


bool Molecule::Impl::isAdjacent(const AtomIndexType a, const AtomIndexType b) const {
  return static_cast<bool>(
    graph::edgeOption(a, b, _adjacencies)
  );
}

bool Molecule::Impl::isSafeToRemoveAtom(const AtomIndexType a) const {
  // A molecule is by definition at least two atoms!
  if(numAtoms() == 2) {
    return false;
  }

  auto removalSafetyData = GraphAlgorithms::getRemovalSafetyData(
    getGraph()
  );

  return removalSafetyData.articulationVertices.count(a) == 0;
}

bool Molecule::Impl::isSafeToRemoveBond(const AtomIndexType a, const AtomIndexType b) const {
  if(auto edgeOption = graph::edgeOption(a, b, _adjacencies)) {
    auto removalSafetyData = GraphAlgorithms::getRemovalSafetyData(
      getGraph()
    );

    return removalSafetyData.bridges.count(
      edgeOption.value()
    ) == 0;
  }

  // Bond does not exist -> it is unsafe to remove
  return false;
}


bool Molecule::Impl::modularCompare(
  const Molecule::Impl& other,
  const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
) const {
  const unsigned thisNumAtoms = numAtoms();

  if(thisNumAtoms != other.numAtoms()) {
    return false;
  }

  auto thisHashes = hashes::generate(getGraph(), getStereocenterList(), comparisonBitmask);
  auto otherHashes = hashes::generate(other.getGraph(), other.getStereocenterList(), comparisonBitmask);

  /* boost isomorphism will allocate a vector of size maxHash, this is dangerous
   * as the maximum hash can be immense, another post-processing step is needed
   * for the calculated hashes to decrease the spatial requirements
   *
   * This maps the hashes to an incremented number
   */
  auto maxHash = hashes::regularize(thisHashes, otherHashes);

  // Where the corresponding index from the other graph is stored
  std::vector<AtomIndexType> indexMap(numAtoms());

  bool isomorphic = boost::isomorphism(
    _adjacencies,
    other._adjacencies,
    boost::make_safe_iterator_property_map(
      indexMap.begin(),
      thisNumAtoms,
      boost::get(boost::vertex_index, _adjacencies)
    ),
    hashes::LookupFunctor(thisHashes),
    hashes::LookupFunctor(otherHashes),
    maxHash,
    boost::get(boost::vertex_index, _adjacencies),
    boost::get(boost::vertex_index, other._adjacencies)
  );

  if(!isomorphic) {
    return false;
  }

  // TODO why are the following two checks still necessary if hashes are used? document why
  if(comparisonBitmask & AtomEnvironmentComponents::BondOrders) {
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
  }

  if(comparisonBitmask & AtomEnvironmentComponents::Symmetries) {
    // Before doing a full equivalence check, peek at the sizes
    if(_stereocenters.size() != other._stereocenters.size()) {
      return false;
    }

    // Check equivalence of the StereocenterLists
    if(
      !_stereocenters.compare(
        other._stereocenters,
        comparisonBitmask
      )
    ) {
      return false;
    }
  }

  return true;
}

unsigned Molecule::Impl::numAtoms() const {
  return boost::num_vertices(_adjacencies);
}

unsigned Molecule::Impl::numBonds() const {
  return boost::num_edges(_adjacencies);
}

RankingInformation Molecule::Impl::rankPriority(
  const AtomIndexType a,
  const std::set<AtomIndexType>& excludeAdjacent,
  const boost::optional<AngstromWrapper>& positionsOption
) const {
  RankingInformation rankingResult;

  // Expects that bond types are set properly, complains otherwise
  rankingResult.ligands = GraphAlgorithms::ligandSiteGroups(
    _adjacencies,
    a,
    excludeAdjacent
  );

  std::string molGraphviz;
#ifndef NDEBUG
  molGraphviz = dumpGraphviz();
#endif

  // Rank the substituents
  auto expandedTree = RankingTree(
    getGraph(),
    getCycleData(),
    getStereocenterList(),
    molGraphviz,
    a,
    excludeAdjacent,
    RankingTree::ExpansionOption::Optimized,
    positionsOption
  );

  rankingResult.sortedSubstituents = expandedTree.getRanked();

  rankingResult.ligandsRanking = RankingInformation::rankLigands(
    rankingResult.ligands,
    rankingResult.sortedSubstituents
  );

  // Find links between them
  rankingResult.links = GraphAlgorithms::substituentLinks(
    _adjacencies,
    getCycleData(),
    a,
    rankingResult.ligands,
    excludeAdjacent
  );

  return rankingResult;
}

std::array<AtomIndexType, 2> Molecule::Impl::vertices(
  const GraphType::edge_descriptor& edge
) const {
  return {
    boost::source(edge, _adjacencies),
    boost::target(edge, _adjacencies)
  };
}

/* Iterators */
RangeForTemporary<GraphType::vertex_iterator> Molecule::Impl::iterateAtoms() const {
  return graph::vertices(_adjacencies);
}

/*! Returns a range-for temporary object allowing c++11 style for loop
 * iteration through an atom's adjacencies
 */
RangeForTemporary<GraphType::adjacency_iterator> Molecule::Impl::iterateAdjacencies(const AtomIndexType a) const {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::iterateAdjacencies: Supplied index is invalid!");
  }

  return graph::adjacents(a, _adjacencies);
}

RangeForTemporary<GraphType::edge_iterator> Molecule::Impl::iterateEdges() const {
  return graph::edges(_adjacencies);
}

RangeForTemporary<GraphType::out_edge_iterator> Molecule::Impl::iterateEdges(const AtomIndexType a) const {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::iterateEdges: Supplied index is invalid!");
  }

  return graph::edges(a, _adjacencies);
}



/* Operators */
RangeForTemporary<GraphType::adjacency_iterator> Molecule::Impl::operator [] (const AtomIndexType a) const {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::operator[]: Supplied index is invalid!");
  }

  return graph::adjacents(a, _adjacencies);
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

Molecule::Molecule(GraphType graph) : _pImpl(
  std::make_unique<Impl>(std::move(graph))
) {}

Molecule::Molecule(
  GraphType graph,
  const AngstromWrapper& positions
) : _pImpl(
  std::make_unique<Impl>(
    std::move(graph),
    positions
  )
) {}

Molecule::Molecule(
  GraphType graph,
  StereocenterList stereocenters
) : _pImpl(
  std::make_unique<Impl>(
    std::move(graph),
    std::move(stereocenters)
  )
) {}

/* Modifiers */
AtomIndexType Molecule::addAtom(
  const Delib::ElementType elementType,
  const AtomIndexType adjacentTo,
  const BondType bondType
) {
  return _pImpl->addAtom(elementType, adjacentTo, bondType);
}

void Molecule::addBond(
  const AtomIndexType a,
  const AtomIndexType b,
  const BondType bondType
) {
  _pImpl->addBond(a, b, bondType);
}

void Molecule::assignStereocenter(
  const AtomIndexType a,
  const boost::optional<unsigned>& assignment
) {
  _pImpl->assignStereocenter(a, assignment);
}

void Molecule::assignStereocenterRandomly(const AtomIndexType a) {
  _pImpl->assignStereocenterRandomly(a);
}

void Molecule::assignStereocenterRandomly(const GraphType::edge_descriptor& e) {
  _pImpl->assignStereocenterRandomly(e);
}

void Molecule::removeAtom(const AtomIndexType a) {
  _pImpl->removeAtom(a);
}

void Molecule::removeBond(
  const AtomIndexType a,
  const AtomIndexType b
) {
  _pImpl->removeBond(a, b);
}

bool Molecule::setBondType(
  const AtomIndexType a,
  const AtomIndexType b,
  const BondType bondType
) {
  return _pImpl->setBondType(a, b, bondType);
}

void Molecule::setElementType(
  const AtomIndexType a,
  const Delib::ElementType elementType
) {
  _pImpl->setElementType(a, elementType);
}

void Molecule::setGeometryAtAtom(
  const AtomIndexType a,
  const Symmetry::Name symmetryName
) {
  _pImpl->setGeometryAtAtom(a, symmetryName);
}


/* Information */
Symmetry::Name Molecule::determineLocalGeometry(
  const AtomIndexType index,
  const RankingInformation& ranking
) const {
  return _pImpl->determineLocalGeometry(index, ranking);
}

std::string Molecule::dumpGraphviz() const {
  return _pImpl->dumpGraphviz();
}

GraphType::edge_descriptor Molecule::edge(
  AtomIndexType a,
  AtomIndexType b
) const {
  return _pImpl->edge(a, b);
}

std::vector<AtomIndexType> Molecule::getAdjacencies(const AtomIndexType a) const {
  return _pImpl->getAdjacencies(a);
}

boost::optional<BondType> Molecule::getBondType(
  const AtomIndexType a,
  const AtomIndexType b
) const {
  return _pImpl->getBondType(a, b);
}

BondType Molecule::getBondType(const GraphType::edge_descriptor& edge) const {
  return _pImpl->getBondType(edge);
}

Cycles Molecule::getCycleData() const {
  return _pImpl->getCycleData();
}

Delib::ElementType Molecule::getElementType(const AtomIndexType index) const {
  return _pImpl->getElementType(index);
}

Delib::ElementTypeCollection Molecule::getElementCollection() const {
  return _pImpl->getElementCollection();
}

const GraphType& Molecule::getGraph() const {
  return _pImpl->getGraph();
}

const StereocenterList& Molecule::getStereocenterList() const {
  return _pImpl->getStereocenterList();
}

unsigned Molecule::getNumAdjacencies(const AtomIndexType a) const {
  return _pImpl->getNumAdjacencies(a);
}

StereocenterList Molecule::inferStereocentersFromPositions(
  const AngstromWrapper& angstromWrapper
) const {
  return _pImpl->inferStereocentersFromPositions(angstromWrapper);
}

bool Molecule::isAdjacent(const AtomIndexType a, const AtomIndexType b) const {
  return _pImpl->isAdjacent(a, b);
}

bool Molecule::isSafeToRemoveAtom(const AtomIndexType a) const {
  return _pImpl->isSafeToRemoveAtom(a);
}

bool Molecule::isSafeToRemoveBond(const AtomIndexType a, const AtomIndexType b) const {
  return _pImpl->isSafeToRemoveBond(a, b);
}

bool Molecule::modularCompare(
  const Molecule& other,
  const temple::Bitmask<AtomEnvironmentComponents>& comparisonBitmask
) const {
  return _pImpl->modularCompare(*other._pImpl, comparisonBitmask);
}

unsigned Molecule::numAtoms() const {
  return _pImpl->numAtoms();
}

unsigned Molecule::numBonds() const {
  return _pImpl->numBonds();
}

RankingInformation Molecule::rankPriority(
  const AtomIndexType a,
  const std::set<AtomIndexType>& excludeAdjacent,
  const boost::optional<AngstromWrapper>& positionsOption
) const {
  return _pImpl->rankPriority(a, excludeAdjacent, positionsOption);
}

std::array<AtomIndexType, 2> Molecule::vertices(
  const GraphType::edge_descriptor& edge
) const {
  return _pImpl->vertices(edge);
}

RangeForTemporary<GraphType::vertex_iterator> Molecule::iterateAtoms() const {
  return _pImpl->iterateAtoms();
}

RangeForTemporary<GraphType::adjacency_iterator> Molecule::iterateAdjacencies(const AtomIndexType a) const {
  return _pImpl->iterateAdjacencies(a);
}

RangeForTemporary<GraphType::edge_iterator> Molecule::iterateEdges() const {
  return _pImpl->iterateEdges();
}

RangeForTemporary<GraphType::out_edge_iterator> Molecule::iterateEdges(const AtomIndexType a) const {
  return _pImpl->iterateEdges(a);
}

/* Operators */
RangeForTemporary<GraphType::adjacency_iterator> Molecule::operator [] (const AtomIndexType a) const {
  return _pImpl->operator[](a);
}

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
  const auto& stereocenters = molecule.getStereocenterList();

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
