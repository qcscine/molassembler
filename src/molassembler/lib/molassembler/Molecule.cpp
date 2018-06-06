#include "temple/Optionals.h"

#include "boost/range/join.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/graph_utility.hpp"

#include "Delib/Constants.h"

#include "chemical_symmetries/ConstexprProperties.h"

#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"

#include "detail/CommonTrig.h"
#include "detail/MolGraphWriter.h"
#include "CNStereocenter.h"
#include "EZStereocenter.h"
#include "GraphAlgorithms.h"
#include "GraphHelpers.h"
#include "LocalGeometryModel.h"
#include "Log.h"
#include "Molecule.h"
#include "RankingTree.h"
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

  /*! Returns if an edge could be an EZStereocenter with multiple assignments
   *
   * Criteria applied are:
   * - Bond type must be double
   * - 2-3 non-eta bonds for each edge vertex
   */
  bool _isEZStereocenterCandidate(const GraphType::edge_descriptor& edgeIndex) const;

  //! Returns whether the specified index is valid or not
  bool _isValidIndex(const AtomIndexType index) const;

  /*!
   * Returns a list of edge indices where each endpoint has 1 or two additional
   * substituent besides the edge neighbor
   */
  std::vector<EdgeIndexType> _getEZStereocenterCandidates() const;

  //! Updates the molecule's StereocenterList after a graph modification
  void _propagateGraphChange();


//!@name Constructors
//!@{
  //! Default constructor
  Impl() noexcept;

  //! Diatomic constructor
  Impl(
    const Delib::ElementType a,
    const Delib::ElementType b,
    const BondType bondType
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
    const Delib::ElementType elementType,
    const AtomIndexType adjacentTo,
    const BondType bondType
  );

  //! Adds a bond between existing atoms.
  void addBond(
    const AtomIndexType a,
    const AtomIndexType b,
    const BondType bondType
  );

  /*! Sets the stereocenter assignment at a particular atom
   *
   * This sets the stereocenter assignment at a specific atom index. For this,
   * a stereocenter must be instantiated and contained in the StereocenterList
   * returned by getStereocenterList(). The supplied assignment must be either
   * boost::none or smaller than stereocenterPtr->numAssignments().
   *
   * @note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenter(
    const AtomIndexType a,
    const boost::optional<unsigned>& assignment
  );

  /*! Assigns a stereocenter stereopermutation at random
   *
   * This sets the stereocetner assignment at a specific index, taking relative
   * statistical occurence weights of each stereopermutation into account. For
   * this, a stereocenter must be instantiated and contained in the
   * StereocenterList returned by getStereocenterList().
   *
   * @note Although molecules in which this occurs are infrequent, consider the
   * StereocenterList you have accessed prior to calling this function and
   * particularly any iterators thereto invalidated. This is because an
   * assignment change can trigger a ranking change, which can in turn lead
   * to the introduction of new stereocenters or the removal of old ones.
   */
  void assignStereocenterRandomly(const AtomIndexType a);

  /*! Removes an atom from the graph, including bonds to it.
   *
   * Removes an atom from the molecular graph, including bonds to the atom,
   * after checking that removing it is safe, i.e. the removal does not
   * disconnect the graph.
   *
   * @throws if the supplied index is invalid or isSafeToRemoveAtom returns false.
   */
  void removeAtom(const AtomIndexType a);

  /*!
   * Removes an atom after checking if removing that bond is safe, i.e. does not
   * disconnect the graph. An example of bonds that can always be removed are
   * ring-closing bonds, since they never disconnect the molecular graph.
   *
   * @throws if isSafeToRemoveBond returns false.
   *
   * @note It is not safe to remove a bond just because one of the involved
   * atoms is terminal, since that atom would then be disconnected from the
   * rest of the molecule. This function merely removes a bond from the graph.
   * It is, however, considered safe to remove the terminal vertex, which
   * involves removing the bond to it.
   */
  void removeBond(const AtomIndexType a, const AtomIndexType b);

  //! Changes an existing bond's type
  bool setBondType(
    const AtomIndexType a,
    const AtomIndexType b,
    const BondType bondType
  );

  //! Changes an existing atom's element type
  void setElementType(
    const AtomIndexType a,
    const Delib::ElementType elementType
  );

  /*! Sets the local geometry at an atom index
   *
   * This sets the local geometry at a specific atom index. There are a number
   * of cases that this function treats differently, besides faulty arguments:
   * If there is already a CNStereocenter instantiated at this atom index, its
   * underlying symmetry is altered. If there is no CNStereocenter at
   * this index, one is instantiated. In all cases, new or modified
   * stereocenters are default-assigned if there is only one possible
   * assignment.
   * @throws if
   *   - the supplied atomic index is invalid
   *   - there is an EZStereocenter at that index
   *   - or the provided symmetry is a different size than that of an existing
   *     CNStereocenter or the expected symmetry
   */
  void setGeometryAtAtom(
    const AtomIndexType a,
    const Symmetry::Name symmetryName
  );
//!@}

//!@name Information
//!@{
  /*! Determines what the local geometry at a non-terminal atom ought to be
   *
   * Returns the expected symmetry name at a non-terminal atom.
   * @throws if the supplied atomic index is invalid
   */
  Symmetry::Name determineLocalGeometry(
    const AtomIndexType index,
    const RankingInformation& ranking
  ) const;

  //! Returns a graphivz string representation of the molecule
  std::string dumpGraphviz() const;

  /*! Fetches the atomic indices of vertices adjacent to a particular index
   *
   * Fetches the atomic indices of vertices adjacent to a particular index.
   * @throws if the supplied atomic index is invalid
   */
  std::vector<AtomIndexType> getAdjacencies(const AtomIndexType a) const;

  //! Fetches the optional bond type between two atom indices
  boost::optional<BondType> getBondType(
    const AtomIndexType a,
    const AtomIndexType b
  ) const;

  BondType getBondType(const GraphType::edge_descriptor& edge) const;

  Cycles getCycleData() const;

  /*! Returns the element type of an atomic index
   *
   * Returns the element type of an atomic index
   * @throws if the atomic index is invalid
   */
  Delib::ElementType getElementType(const AtomIndexType index) const;

  //! Returns a collection detailing all element types
  Delib::ElementTypeCollection getElementCollection() const;

  //! Provides read-only access to the graph member
  const GraphType& getGraph() const;

  //! Provides read-only access to the list of stereocenters
  const StereocenterList& getStereocenterList() const;

  //! Returns the number of adjacencies of an atomic position
  unsigned getNumAdjacencies(const AtomIndexType a) const;

  StereocenterList inferStereocentersFromPositions(const AngstromWrapper& angstromWrapper) const;

  //! Checks if two atomic indices are connected by a bond
  bool isAdjacent(
    const AtomIndexType a,
    const AtomIndexType b
  ) const;

  //! An atom is considered removable if it isn't an articulation vertex
  bool isSafeToRemoveAtom(const AtomIndexType a) const;

  //! A bond is considered removable if it isn't a bridge edge
  bool isSafeToRemoveBond(const AtomIndexType a, const AtomIndexType b) const;

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
    const AtomIndexType a,
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
    const AtomIndexType a
  ) const;

  /*! Returns a range-for temporary object allowing c++11-style for loop
   * iteration through edges
   */
  RangeForTemporary<GraphType::edge_iterator> iterateEdges() const;

  /*! Returns a range-for temporary object allowing c++11-style for loop
   * iteration through edges around a specific atom
   */
  RangeForTemporary<GraphType::out_edge_iterator> iterateEdges(const AtomIndexType a) const;
//!@}


//!@name Operators
//!@{
  //! Returns the adjacencies of the specified atom index
  RangeForTemporary<GraphType::adjacency_iterator> operator [] (
    const AtomIndexType a
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

  /* TODO
   * - Will need refinement to not instantiate EZStereocenters in small cycles
   *   (up to a preset size, maybe around 8 or so?)
   */
  // Find EZStereocenters
  for(const auto& edgeIndex : _getEZStereocenterCandidates()) {
    AtomIndexType source = boost::source(edgeIndex, _adjacencies),
                  target = boost::target(edgeIndex, _adjacencies);

    // Construct a Stereocenter here
    auto newStereocenter = std::make_shared<Stereocenters::EZStereocenter>(
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

  Cycles cycleData = getCycleData();
  // Find CNStereocenters
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
    auto newStereocenter = std::make_shared<Stereocenters::CNStereocenter>(
      _adjacencies,
      determineLocalGeometry(candidateIndex, localRanking),
      candidateIndex,
      localRanking
    );

    if(newStereocenter -> numAssignments() == 1) {
      newStereocenter -> assign(0);
    }

    if(
      !disregardStereocenter(
        *newStereocenter,
        getElementType(candidateIndex),
        cycleData,
        Options::temperatureRegime
      )
    ) {
      stereocenterList.add(
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

bool Molecule::Impl::_isEZStereocenterCandidate(const GraphType::edge_descriptor& edgeIndex) const {
  auto numNonEtaAdjacencies = [&](const AtomIndexType a) -> unsigned {
    unsigned nonEta = 0;

    for(const auto& edgeIndex : iterateEdges(a)) {
      if(_adjacencies[edgeIndex].bondType != BondType::Eta) {
        nonEta += 1;
      }
    }

    return nonEta;
  };

  auto source = boost::source(edgeIndex, _adjacencies),
       target = boost::target(edgeIndex, _adjacencies);

  auto sourceAdjacencies = numNonEtaAdjacencies(source),
       targetAdjacencies = numNonEtaAdjacencies(target);

  return (
    _adjacencies[edgeIndex].bondType == BondType::Double
    && 2 <= sourceAdjacencies
    && sourceAdjacencies <= 3
    && 2 <= targetAdjacencies
    && targetAdjacencies <= 3
  );
}

std::vector<EdgeIndexType> Molecule::Impl::_getEZStereocenterCandidates() const {
  std::vector<EdgeIndexType> candidates;

  for(const auto& edgeIndex : iterateEdges()) {
    if(_isEZStereocenterCandidate(edgeIndex)) {
      candidates.push_back(edgeIndex);
    }
  }

  return candidates;
}

void Molecule::Impl::_propagateGraphChange() {
  /* Two cases: If the StereocenterList is empty, we can just use detect to
   * find any new stereocenters in the Molecule.
   *
   * In the other case, we have to recheck absolutely everywhere. If ranking
   * was affected and the stereocenter has a set assignment, we need to find
   * the assignment that the previous ranking represented spatially in the new
   * set of assignments and assign the stereocenter to that.
   */
  if(_stereocenters.empty()) {
    _stereocenters = _detectStereocenters();
  } else {
    GraphAlgorithms::findAndSetEtaBonds(_adjacencies);

    // EZStereocenters first
    for(const auto& edgeIndex : iterateEdges()) {
      auto source = boost::source(edgeIndex, _adjacencies),
           target = boost::target(edgeIndex, _adjacencies);

      // Is there already a stereocenter on both edge vertices?
      if(_stereocenters.involving(source) && _stereocenters.involving(target)) {
        // Is it the same, and an EZStereocenter?
        if(
          _stereocenters.at(source) == _stereocenters.at(target)
          && _stereocenters.at(source)->type() == Stereocenters::Type::EZStereocenter
        ) {
          // Is it even possible for it to be a candidate anymore?
          if(_isEZStereocenterCandidate(edgeIndex)) {
            // Re-rank, and adapt it to a new ranking
            std::dynamic_pointer_cast<Stereocenters::EZStereocenter>(
              _stereocenters.at(source)
            ) -> propagateGraphChange(
              rankPriority(source, {target}),
              rankPriority(target, {source})
            );

            // If this EZStereocenter has only one assignment now, remove it
            if(_stereocenters.at(source) -> numAssignments() == 1) {
              _stereocenters.remove(source);
            }
          } else {
            // It cannot be an EZStereocenter anymore, so drop it from the list
            _stereocenters.remove(source);
          }
        }
        /* If the stereocenters are NOT the same, in which case both would
         * have to be CNStereocenters, then this edge cannot be an
         * EZStereocenter. Do nothing.
         */
      } else if(
        !_stereocenters.involving(source)
        && !_stereocenters.involving(target)
        && _isEZStereocenterCandidate(edgeIndex)
      ) {
        auto newStereocenterPtr = std::make_shared<Stereocenters::EZStereocenter>(
          source,
          rankPriority(source, {target}),
          target,
          rankPriority(target, {source})
        );

        if(newStereocenterPtr -> numAssignments() == 2) {
          _stereocenters.add(
            std::move(newStereocenterPtr)
          );
        }
      }
      /* If there is a stereocenter on one vertex, but not the other, then it
       * should be a CNStereocenter, in which case we best do nothing.
       */
    }

    Cycles cycleData = getCycleData();
    // Now CNStereocenters
    for(
      AtomIndexType candidateIndex = 0;
      candidateIndex < numAtoms();
      ++candidateIndex
    ) {
      if(_stereocenters.involving(candidateIndex)) {
        if(_stereocenters.at(candidateIndex)->type() == Stereocenters::Type::CNStereocenter) {
          // Is it possible for this atom to continue to be a CNStereocenter?
          auto CNStereocenterPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
            _stereocenters.at(candidateIndex)
          );

          RankingInformation localRanking = rankPriority(candidateIndex);

          // If the atom is now terminal, remove the stereocenter
          if(localRanking.ligands.size() <= 1) {
            _stereocenters.remove(candidateIndex);
            continue;
          }

          // Propagate the state of the stereocenter to the new ranking
          CNStereocenterPtr -> propagateGraphChange(
            _adjacencies,
            localRanking
          );

          /* If the modified stereocenter has only one assignment and is
           * unassigned due to the graph change, default-assign it
           */
          if(CNStereocenterPtr -> numAssignments() == 1) {
            if(CNStereocenterPtr -> assigned() == boost::none) {
              CNStereocenterPtr -> assign(0);
            }
          }

          if(
            disregardStereocenter(
              *CNStereocenterPtr,
              getElementType(candidateIndex),
              cycleData,
              Options::temperatureRegime
            )
          ) {
            // Since this index cannot be a CNStereocenter anymore, drop it
            _stereocenters.remove(candidateIndex);
          }
        }
        /* If the stereocenter at that candidate index is an EZStereocenter, do
         * nothing
         */
      } else {
        // No stereocenter yet
        auto localRanking = rankPriority(candidateIndex);

        // Skip terminal atoms
        if(localRanking.ligands.size() <= 1) {
          continue;
        }

        auto newStereocenterPtr = std::make_shared<Stereocenters::CNStereocenter>(
          _adjacencies,
          determineLocalGeometry(candidateIndex, localRanking),
          candidateIndex,
          localRanking
        );

        if(newStereocenterPtr -> numAssignments() == 1) {
          newStereocenterPtr -> assign(0);
        }

        if(
          !disregardStereocenter(
            *newStereocenterPtr,
            getElementType(candidateIndex),
            cycleData,
            Options::temperatureRegime
          )
        ) {
          _stereocenters.add(
            std::move(newStereocenterPtr)
          );
        }
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

Molecule::Impl::Impl(GraphType graph)
: _adjacencies(std::move(graph))
{
  GraphAlgorithms::findAndSetEtaBonds(_adjacencies);
  _stereocenters = _detectStereocenters();
  _ensureModelInvariants();
}

Molecule::Impl::Impl(
  GraphType graph,
  const AngstromWrapper& angstromWrapper
) : _adjacencies(std::move(graph))
{
  GraphAlgorithms::findAndSetEtaBonds(_adjacencies);
  _stereocenters = inferStereocentersFromPositions(angstromWrapper);
  _ensureModelInvariants();
}

Molecule::Impl::Impl(
  GraphType graph,
  StereocenterList stereocenters
) : _adjacencies(std::move(graph)),
    _stereocenters(std::move(stereocenters))
{
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
    if(_stereocenters.involving(toIndex)) {
      if(_stereocenters.at(toIndex)->type() == Stereocenters::Type::CNStereocenter) {
        auto localRanking = rankPriority(toIndex);

        std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
          _stereocenters.at(toIndex)
        ) -> addSubstituent(
          _adjacencies,
          addedIndex,
          localRanking,
          determineLocalGeometry(toIndex, localRanking),
          Options::chiralStatePreservation
        );
      } else {
        // Adding this new adjacency invalidates the EZStereocenter there
        if(getNumAdjacencies(toIndex) > 2) {
          _stereocenters.remove(toIndex);
        } else {
          auto EZPtr = std::dynamic_pointer_cast<Stereocenters::EZStereocenter>(
            _stereocenters.at(toIndex)
          );

          // What is the other central index?
          AtomIndexType otherCenter = (
            EZPtr -> involvedAtoms().at(0) == toIndex
            ? EZPtr -> involvedAtoms().at(1)
            : EZPtr -> involvedAtoms().at(0)
          );

          EZPtr -> addSubstituent(
            toIndex,
            rankPriority(toIndex, {otherCenter})
          );
        }
      }
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

  if(!_stereocenters.involving(a)) {
    throw std::logic_error("assignStereocenter: No stereocenter at this index!");
  }

  auto stereocenterPtr = _stereocenters.at(a);

  if(assignment >= stereocenterPtr -> numAssignments()) {
    throw std::logic_error("assignStereocenter: Invalid assignment index!");
  }

  stereocenterPtr -> assign(assignment);

  // A reassignment can change ranking! See the RankingTree tests
  _propagateGraphChange();
}

void Molecule::Impl::assignStereocenterRandomly(const AtomIndexType a) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::assignStereocenterRandomly: Supplied index is invalid!");
  }

  if(!_stereocenters.involving(a)) {
    throw std::logic_error("assignStereocenterRandomly: No stereocenter at this index!");
  }

  _stereocenters.at(a) -> assignRandom();

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
  if(_stereocenters.involving(a)) {
    _stereocenters.remove(a);
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
    if(_stereocenters.involving(indexToUpdate)) {
      if(_stereocenters.at(indexToUpdate) -> type() == Stereocenters::Type::CNStereocenter) {
        auto localRanking = rankPriority(indexToUpdate);

        /* If the index on which the CNStereocenter is placed becomes terminal,
         * drop the stereocenter
         */
        if(localRanking.ligands.size() <= 1) {
          _stereocenters.remove(indexToUpdate);
          continue;
        }

        std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
          _stereocenters.at(indexToUpdate)
        ) -> removeSubstituent(
          _adjacencies,
          Stereocenters::Stereocenter::removalPlaceholder,
          localRanking,
          determineLocalGeometry(indexToUpdate, localRanking),
          Options::chiralStatePreservation
        );
      } else {
        std::dynamic_pointer_cast<Stereocenters::EZStereocenter>(
          _stereocenters.at(indexToUpdate)
        ) -> removeSubstituent(
          indexToUpdate,
          Stereocenters::Stereocenter::removalPlaceholder
        );
      }
    }
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
  if(edgePair.second) {
    boost::remove_edge(edgePair.first, _adjacencies);

    /* If there is an EZStereocenter on this edge, we have to drop it explicitly,
     * since _propagateGraphChange cannot iterate over a now-removed edge.
     */
    if(
      _stereocenters.involving(a)
      && _stereocenters.involving(b)
      && _stereocenters.at(a) == _stereocenters.at(b)
    ) {
      _stereocenters.remove(a);
    }

    // Notify all immediately adjacent stereocenters of the removal
    auto notifyRemoval = [this](
      const auto& indexToUpdate,
      const auto& removedIndex
    ) {
      if(_stereocenters.involving(indexToUpdate)) {
        if(_stereocenters.at(indexToUpdate) -> type() == Stereocenters::Type::CNStereocenter) {
          auto localRanking = rankPriority(indexToUpdate);

          // In case the CNS central atom becomes terminal, just drop the stereocenter
          if(localRanking.ligands.size() <= 1) {
            _stereocenters.remove(indexToUpdate);
            return;
          }

          std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
            _stereocenters.at(indexToUpdate)
          ) -> removeSubstituent(
            _adjacencies,
            removedIndex,
            localRanking,
            determineLocalGeometry(indexToUpdate, localRanking),
            Options::chiralStatePreservation
          );
        } else {
          std::dynamic_pointer_cast<Stereocenters::EZStereocenter>(
            _stereocenters.at(indexToUpdate)
          ) -> removeSubstituent(
            indexToUpdate,
            removedIndex
          );
        }
      }
    };

    notifyRemoval(a, b);
    notifyRemoval(b, a);

    /* All other cases, where there may be EZStereocenters or CNStereocenters
     * on a or b, should be handled correctly by _propagateGraphChange.
     */

    _propagateGraphChange();
  } else {
    throw std::logic_error("That bond does not exist!");
  }
}

bool Molecule::Impl::setBondType(
  const AtomIndexType a,
  const AtomIndexType b,
  const BondType bondType
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

  if(_stereocenters.involving(a)) {
    if(_stereocenters.at(a)->type() == Stereocenters::Type::CNStereocenter) {
      auto CNSPointer = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
        _stereocenters.at(a)
      );

      if(
        Symmetry::size(CNSPointer->getSymmetry())
        == Symmetry::size(symmetryName)
      ) {
        CNSPointer->setSymmetry(symmetryName, _adjacencies);
        _propagateGraphChange();
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
    RankingInformation localRanking = rankPriority(a);
    const Symmetry::Name expectedSymmetry = determineLocalGeometry(a, localRanking);

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
        _adjacencies,
        symmetryName,
        a,
        localRanking
      );

      // Default-assign stereocenters with only one assignment
      if(newStereocenterPtr->numAssignments() == 1) {
        newStereocenterPtr->assign(0u);
      }

      _stereocenters.add(
        std::move(newStereocenterPtr)
      );

      _propagateGraphChange();
    }
  }
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

  for(const auto& edgeIndex : _getEZStereocenterCandidates()) {
    auto source = graph::source(edgeIndex, _adjacencies),
         target = graph::target(edgeIndex, _adjacencies);

    // Construct a Stereocenter here
    auto newStereocenter = std::make_shared<Stereocenters::EZStereocenter>(
      source,
      rankPriority(source, {target}, angstromWrapper),
      target,
      rankPriority(target, {source}, angstromWrapper)
    );

    newStereocenter -> fit(angstromWrapper);

    if(newStereocenter -> numAssignments() == 2) {
      stereocenters.add(
        std::move(newStereocenter)
      );
    }
  }

  Cycles cycleData = getCycleData();
  /* Add a CNStereocenter everywhere where the symmetry yielding the best fit is
   * not the one that Molecule's determineLocalGeometry gets and where we
   * can fully determine a Stereocenter's assignment from the positions
   */
  for(unsigned candidateIndex = 0; candidateIndex < numAtoms(); candidateIndex++) {
    // Skip unsuitable atoms and ones that already have a stereocenter
    if(stereocenters.involving(candidateIndex)) {
      continue;
    }

    RankingInformation localRanking = rankPriority(candidateIndex, {}, angstromWrapper);

    // Skip terminal atoms
    if(localRanking.ligands.size() <= 1) {
      continue;
    }

    const Symmetry::Name expectedGeometry = determineLocalGeometry(candidateIndex, localRanking);

    // Construct it
    auto stereocenterPtr = std::make_shared<Stereocenters::CNStereocenter>(
      _adjacencies,
      expectedGeometry,
      candidateIndex,
      localRanking
    );

    pickyFit(
      *stereocenterPtr,
      _adjacencies,
      angstromWrapper,
      expectedGeometry
    );

    if(
      !disregardStereocenter(
        *stereocenterPtr,
        graph::elementType(candidateIndex, _adjacencies),
        cycleData,
        Options::temperatureRegime
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
  using HashType = hashes::AtomEnvironmentHashType;
  using ComparisonComponents = AtomEnvironmentComponents;

  const unsigned thisNumAtoms = numAtoms();

  if(thisNumAtoms != other.numAtoms()) {
    return false;
  }

  auto generateHashes = [&comparisonBitmask](const Molecule::Impl& mol) -> std::vector<HashType> {

    std::vector<HashType> hashes;

    const AtomIndexType N = mol.numAtoms();
    hashes.reserve(N);

    std::vector<BondType> bonds;
    bonds.reserve(Symmetry::constexprProperties::maxSymmetrySize);

    for(AtomIndexType i = 0; i < N; ++i) {
      if(comparisonBitmask & ComparisonComponents::BondOrders) {
        bonds.clear();

        for(const auto& edge : mol.iterateEdges(i)) {
          bonds.emplace_back(mol.getGraph()[edge].bondType);
        }

        std::sort(
          bonds.begin(),
          bonds.end()
        );
      }

      boost::optional<Symmetry::Name> symmetryNameOption;
      boost::optional<unsigned> assignmentOption;

      if(mol.getStereocenterList().involving(i)) {
        const auto& stereocenterPtr = mol.getStereocenterList().at(i);
        if(stereocenterPtr->type() == Stereocenters::Type::EZStereocenter) {
          symmetryNameOption = Symmetry::Name::TrigonalPlanar;
        } else {
          symmetryNameOption = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
            stereocenterPtr
          ) -> getSymmetry();
        }

        assignmentOption = stereocenterPtr -> assigned();
      }

      hashes.emplace_back(
        hashes::atomEnvironment(
          comparisonBitmask,
          mol.getElementType(i),
          bonds,
          symmetryNameOption,
          assignmentOption
        )
      );
    }

    return hashes;
  };

  auto thisHashes = generateHashes(*this);
  auto otherHashes = generateHashes(other);

  /* boost isomorphism will allocate a vector of size maxHash, this is dangerous
   * as the maximum hash can be immense, another post-processing step is needed
   * for the calculated hashes to decrease the spatial requirements
   */
  std::unordered_map<HashType, HashType> reductionMapping;
  HashType counter = 0;

  // Generate mapping from hash values to integer-incremented reduction
  for(const auto& hash : boost::range::join(thisHashes, otherHashes)) {
    if(reductionMapping.count(hash) == 0) {
      reductionMapping.emplace(
        hash,
        counter
      );

      ++counter;
    }
  }

  // Reduce the hash representations
  for(auto& hash : boost::range::join(thisHashes, otherHashes)) {
    hash = reductionMapping.at(hash);
  }

  auto maxHash = counter;

  if(maxHash > std::numeric_limits<GraphType::vertices_size_type>::max()) {
    throw std::logic_error(
      "Number of distinct atom environment hashes exceeds limits of boost "
      " graph's isomorphism algorithm type used to store it!"
    );
  }

  // This explicit form is needed instead of a lambda for boost's concept checks
  struct HashLookup {
    const std::vector<HashType>* const hashes;

    using argument_type = AtomIndexType;
    using result_type = HashType;

    HashLookup(const std::vector<HashType>& hashes) : hashes(&hashes) {}

    HashType operator() (const AtomIndexType i) const {
      return hashes->at(i);
    }
  };

  /* Ensure this matches the required boost concept. Although the documentation
   * asks for the UnaryFunction concept, in reality it requires
   * AdaptableUnaryFunction with typedefs for argument and result.
   */
  BOOST_CONCEPT_ASSERT((
    boost::AdaptableUnaryFunctionConcept<HashLookup, HashType, AtomIndexType>
  ));

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
    HashLookup(thisHashes),
    HashLookup(otherHashes),
    maxHash,
    boost::get(boost::vertex_index, _adjacencies),
    boost::get(boost::vertex_index, other._adjacencies)
  );

  if(!isomorphic) {
    return false;
  }

  if(comparisonBitmask & ComparisonComponents::BondOrders) {
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

  if(comparisonBitmask & ComparisonComponents::Symmetries) {
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
        ) {
          return false;
        }

        if(
          (comparisonBitmask & ComparisonComponents::Stereopermutations)
          && thisCNSPtr->assigned() != otherCNSPtr->assigned()
        ) {
          return false;
        }
      } else {
        const auto otherCentralAtoms = temple::map(
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
        if(stereocenterPtr->numAssignments() != otherPtr->numAssignments()) {
          return false;
        }

        if(
          (comparisonBitmask & ComparisonComponents::Stereopermutations)
          && stereocenterPtr->assigned() != otherPtr->assigned()
        ) {
          return false;
        }
      }
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
  using ComparisonComponents = AtomEnvironmentComponents;
  // Operator == performs the most strict comparison possible
  return modularCompare(
    other,
    temple::make_bitmask(ComparisonComponents::ElementTypes)
      | ComparisonComponents::BondOrders
      | ComparisonComponents::Symmetries
      | ComparisonComponents::Stereopermutations
  );
}

bool Molecule::Impl::operator != (const Impl& other) const {
  return !(*this == other);
}

/* Molecule interface to Impl call forwards */
Molecule::Molecule() noexcept : _pImpl(
  std::make_unique<Impl>()
) {}

Molecule::Molecule(Molecule&& other) = default;
Molecule& Molecule::operator = (Molecule&& rhs) = default;

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
  const AngstromWrapper& angstromWrapper
) : _pImpl(
  std::make_unique<Impl>(
    std::move(graph),
    angstromWrapper
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
  if(!molecule.getStereocenterList().empty()) {
    os << "Stereocenter information:\n";

    for(const auto& stereocenterPtr: molecule.getStereocenterList()) {
      os << stereocenterPtr -> info() << std::endl;
    }
  }

  return os;
}
