#include "temple/Optionals.h"

#include "boost/range/join.hpp"
#include "boost/graph/biconnected_components.hpp"
#include "boost/graph/graphviz.hpp"
#include "boost/graph/isomorphism.hpp"
#include "boost/graph/graph_utility.hpp"

#include "chemical_symmetries/ConstexprProperties.h"

#include "temple/Containers.h"
#include "temple/constexpr/Numeric.h"

#include "CNStereocenter.h"
#include "CommonTrig.h"
#include "EZStereocenter.h"
#include "GraphAlgorithms.h"
#include "Log.h"
#include "Molecule.h"
#include "MolGraphWriter.h"
#include "RankingTree.h"

namespace molassembler {

/* Molecule implementation ---------------------------------------------------*/
/* "Global" options */
TemperatureRegime Molecule::temperatureRegime = TemperatureRegime::High;
ChiralStatePreservation Molecule::chiralStatePreservation = ChiralStatePreservation::EffortlessAndUnique;

/* Static functions */
Delib::BondOrderCollection Molecule::uffBondOrders(
  const Delib::AtomCollection& atomCollection
) {
  const int N = atomCollection.size();

  auto bondOrders = Delib::BondOrderCollection::createEmpty(N);

  for(int i = 0; i < N; ++i) {
    for(int j = i + 1; j < N; ++j) {
      bondOrders.setOrder(
        i,
        j,
        Bond::calculateBondOrder(
          atomCollection.getElement(i),
          atomCollection.getElement(j),
          (
            atomCollection.getPosition(j)
            - atomCollection.getPosition(i)
          ).norm()
        )
      );
    }
  }

  return bondOrders;
}

Molecule::PseudoHashType Molecule::hashAtomEnvironment(
  const Delib::ElementType& elementType,
  const std::vector<BondType>& sortedBonds,
  boost::optional<Symmetry::Name> symmetryNameOptional,
  boost::optional<unsigned> assignedOptional
) {
  static_assert(
    (
      7
      + 4 * Symmetry::constexprProperties::maxSymmetrySize
      + temple::Math::ceil(
        temple::Math::log(Symmetry::nSymmetries + 1.0, 2.0)
      )
    ) < 64,
    "Element type, bond and symmetry information no longer fit into a 64-bit unsigned integer"
  );


  /* First 8 bits of the 64 bit unsigned number are from the element type
   *
   * Biggest is Cn, which has a value of 112 -> Fits in 7 bits (2^7 = 128)
   */
  PseudoHashType value = static_cast<PseudoHashType>(elementType);

  /* Bonds have 8 possible values currently, plus None is 9
   * -> fits into 4 bits (2^4 = 16)
   *
   * (You could make an argument for removing Sextuple from the list of bond
   * types and fitting this precisely into 3 bits only)
   *
   * So, left shift by 7 bits (so there are 7 zeros on the right in the bit
   * representation) plus the current bond number multiplied by 4 to place a
   * maximum of 8 bonds (maximum symmetry size currently)
   *
   * This occupies 4 * 8 = 32 bits.
   */
  unsigned bondNumber = 0;
  for(const auto& bond : sortedBonds) {
    value += (
      // No bond is represented by 0, while the remaining bond types are shifted
      static_cast<PseudoHashType>(bond) + 1
    ) << (7 + 4 * bondNumber);

    ++bondNumber;
  }

  if(symmetryNameOptional) {
    /* We add symmetry information on non-terminal atoms. There are currently
     * 16 symmetries, plus None is 17, which fits into 5 bits (2^5 = 32)
     */
    value += (static_cast<PseudoHashType>(symmetryNameOptional.value()) + 1) << 39;

    /* The remaining space 64 - (7 + 32 + 5) = 20 bits is used for the current
     * permutation. In that space, we can store up to 2^20 - 2 > 1e5
     * permutations (one (0) for no stereocenter, one (1) for unassigned),
     * which ought to be plenty of bit space. Maximally asymmetric square
     * antiprismatic has around 6k permutations, which fits into 13 bits.
     */
    PseudoHashType permutationValue;
    if(assignedOptional) {
      permutationValue = assignedOptional.value() + 2;
    } else {
      permutationValue = 1;
    }
    value += static_cast<PseudoHashType>(permutationValue) << 44;
  }

  return value;
}


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
    auto newStereocenter = std::make_shared<Stereocenters::EZStereocenter>(
      source,
      rankPriority(source, {target}),
      target,
      rankPriority(target, {source})
    );

    if(newStereocenter -> numStereopermutations() == 2) {
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
        *this,
        determineLocalGeometry(candidateIndex),
        candidateIndex,
        rankPriority(candidateIndex)
      );

      if(newStereocenter -> numStereopermutations() > 1) {
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

bool Molecule::_isCNStereocenterCandidate(const AtomIndexType& atomIndex) const {
  auto numAdjacencies = getNumAdjacencies(atomIndex);

  if(numAdjacencies < 3) {
    return false;
  }

  if(temperatureRegime == TemperatureRegime::High) {
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

bool Molecule::_isEZStereocenterCandidate(const GraphType::edge_descriptor& edgeIndex) const {
  auto numNonEtaAdjacencies = [&](const AtomIndexType& a) -> unsigned {
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

std::vector<EdgeIndexType> Molecule::_getEZStereocenterCandidates() const {
  std::vector<EdgeIndexType> candidates;

  for(
    const auto& edgeIndex :
    RangeForTemporary<GraphType::edge_iterator>(
      boost::edges(_adjacencies)
    )
  ) {
    if(_isEZStereocenterCandidate(edgeIndex)) {
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
  AtomIndexType centralAtom = stereocenter.involvedAtoms().front();

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
      *this,
      positions,
      {Symmetry::Name::Seesaw}
    );
  } else {
    stereocenter.fit(*this, positions);
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

void Molecule::_propagateGraphChange() {
  /* Two cases: If the StereocenterList is empty, we can just use detect to
   * find stereocenters in the Molecule, if there are any new ones.
   *
   * In the other case, we have to recheck everywhere. If ranking was affected
   * and the stereocenter has a set assignment, we need to find the assignment
   * that the previous ranking represented spatially in the new set of
   * assignments and assign the stereocenter to that.
   */
  if(_stereocenters.empty()) {
    _stereocenters = _detectStereocenters();
  } else {
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
            if(_stereocenters.at(source) -> numStereopermutations() == 1) {
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

        if(newStereocenterPtr -> numStereopermutations() == 2) {
          _stereocenters.add(
            std::move(newStereocenterPtr)
          );
        }
      }
      /* If there is a stereocenter on one vertex, but not the other, then it
       * should be a CNStereocenter, in which case we best do nothing.
       */
    }

    // Now CNStereocenters
    for(
      AtomIndexType candidateIndex = 0;
      candidateIndex < numAtoms();
      ++candidateIndex
    ) {
      if(_stereocenters.involving(candidateIndex)) {
        if(_stereocenters.at(candidateIndex)->type() == Stereocenters::Type::CNStereocenter) {
          // Is it possible for this atom to continue to be a CNStereocenter?
          if(_isCNStereocenterCandidate(candidateIndex)) {
            auto CNStereocenterPtr = std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
              _stereocenters.at(candidateIndex)
            );

            // Propagate the state of the stereocenter to the new ranking
            CNStereocenterPtr -> propagateGraphChange(
              *this,
              rankPriority(candidateIndex)
            );

            /* If the modified stereocenter has only one assignment and the
             * determined symmetry adds nothing, remove it
             */
            if(CNStereocenterPtr -> numStereopermutations() == 1) {
              if(CNStereocenterPtr -> getSymmetry() == determineLocalGeometry(candidateIndex)) {
                _stereocenters.remove(candidateIndex);
              }
            }
          } else {
            // Since this index cannot be a CNStereocenter anymore, drop it
            _stereocenters.remove(candidateIndex);
          }
        }
        /* If the stereocenter at that candidate index is an EZStereocenter, do
         * nothing
         */
      } else {
        // No stereocenter yet
        if(_isCNStereocenterCandidate(candidateIndex)) {
          auto newStereocenterPtr = std::make_shared<Stereocenters::CNStereocenter>(
            *this,
            determineLocalGeometry(candidateIndex),
            candidateIndex,
            rankPriority(candidateIndex)
          );

          if(newStereocenterPtr -> numStereopermutations() > 1) {
            _stereocenters.add(
              std::move(newStereocenterPtr)
            );
          }
        }
      }
    }
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

Molecule::Molecule(
  const Delib::AtomCollection& atomCollection,
  const Delib::BondOrderCollection& bondOrders
) : _adjacencies(atomCollection.size()) {
  // Discretize bond orders
  const int N = atomCollection.size();

  for(int i = 0; i < N; ++i) {
    for(int j = i + 1; j < N; ++j) {
      double bondOrder = bondOrders.getOrder(i, j);

      if(bondOrder > 0.5) {
        BondType bond = static_cast<BondType>(
          std::round(bondOrder) - 1
        );

        if(bondOrder > 6.5) {
          bond = BondType::Sextuple;
        }

        auto edgeAddPair = boost::add_edge(i, j, _adjacencies);

        _adjacencies[edgeAddPair.first].bondType = bond;
      }
    }
  }

  _stereocenters = inferStereocentersFromPositions(
    atomCollection.getPositions()
  );
}

Molecule::Molecule(const Delib::AtomCollection& atomCollection)
  : Molecule {atomCollection, uffBondOrders(atomCollection)} {}

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
  // addBond handles the stereocenter update on adjacentTo

  _propagateGraphChange();

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

  if(!edgeAddPair.second) {
    throw std::logic_error("Molecule::addbond: Cannot add a bond where one already is present!");
  }

  _adjacencies[edgeAddPair.first].bondType = bondType;

  auto notifySubstituentAddition = [this](
    const AtomIndexType& toIndex,
    const AtomIndexType& addedIndex
  ) {
    if(_stereocenters.involving(toIndex)) {
      if(_stereocenters.at(toIndex)->type() == Stereocenters::Type::CNStereocenter) {
        std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
          _stereocenters.at(toIndex)
        ) -> addSubstituent(
          *this,
          addedIndex,
          rankPriority(toIndex),
          determineLocalGeometry(toIndex),
          chiralStatePreservation
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

void Molecule::assignStereocenterAtAtom(
  const AtomIndexType& a,
  const boost::optional<unsigned>& assignment
) {
  if(!_isValidIndex(a)) {
    throw std::out_of_range("Molecule::assignStereocenterAtAtom: Supplied index is invalid!");
  }

  if(_stereocenters.involving(a)) {
    auto stereocenterPtr = _stereocenters.at(a);

    if(assignment < stereocenterPtr -> numStereopermutations()) {
      stereocenterPtr -> assign(assignment);

      // A reassignment can change ranking! See the RankingTree tests
      _propagateGraphChange();
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
   * liberally in the stereocenter classes. We have to correct of all of those
   * to ensure that _propagateGraphChange works properly.
   */
  _stereocenters.propagateVertexRemoval(a);

  /* call removeSubstituent on all adjacent stereocenters, with
   * std::numeric_limits<AtomIndexType>::max() as the 'which' parameter,
   * which is what propagateVertexRemoval replaces the removed index with in the
   * stereocenters' internal state
   */
  for(const auto& indexToUpdate : previouslyAdjacentVertices) {
    if(_stereocenters.involving(indexToUpdate)) {
      if(_stereocenters.at(indexToUpdate) -> type() == Stereocenters::Type::CNStereocenter) {
        std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
          _stereocenters.at(indexToUpdate)
        ) -> removeSubstituent(
          *this,
          std::numeric_limits<AtomIndexType>::max(),
          rankPriority(indexToUpdate),
          determineLocalGeometry(indexToUpdate),
          chiralStatePreservation
        );
      } else {
        std::dynamic_pointer_cast<Stereocenters::EZStereocenter>(
          _stereocenters.at(indexToUpdate)
        ) -> removeSubstituent(
          indexToUpdate,
          std::numeric_limits<AtomIndexType>::max()
        );
      }
    }
  }

  _propagateGraphChange();
}

void Molecule::removeBond(
  const AtomIndexType& a,
  const AtomIndexType& b
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
          std::dynamic_pointer_cast<Stereocenters::CNStereocenter>(
            _stereocenters.at(indexToUpdate)
          ) -> removeSubstituent(
            *this,
            removedIndex,
            rankPriority(indexToUpdate),
            determineLocalGeometry(indexToUpdate),
            chiralStatePreservation
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

  _propagateGraphChange();

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
  _propagateGraphChange();
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
        CNSPointer->setSymmetry(*this, symmetryName);
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
        *this,
        symmetryName,
        a,
        rankPriority(a)
      );

      // Default-assign stereocenters with only one assignment
      if(newStereocenterPtr->numStereopermutations() == 1) {
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
Symmetry::Name Molecule::determineLocalGeometry(
  const AtomIndexType& index
) const {
  if(!_isValidIndex(index)) {
    throw std::out_of_range("Molecule::determineLocalGeometry: Supplied index is invalid!");
  }

  if(getNumAdjacencies(index) <= 1) {
    throw std::logic_error(
      "Molecule::determineLocalGeometry: No geometries exist for terminal atoms"
    );
  }

  auto ligandsVector = _reduceToLigandTypes(index);

  // TODO this below is invalid for metals!
  unsigned nSites = getNumAdjacencies(index);
  int formalCharge = 0;

  auto symmetryOptional = LocalGeometry::vsepr(
    getElementType(index),
    nSites,
    ligandsVector,
    formalCharge
  ) | temple::callIfNone(LocalGeometry::firstOfSize, nSites);

  if(!symmetryOptional) {
    throw std::logic_error(
      "Could not determine a geometry! Perhaps you have more substituents "
      "than the largest symmety can handle?"
    );
  }

  return symmetryOptional.value();
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

Delib::ElementTypeCollection Molecule::getElementCollection() const {
  AtomIndexType N = numAtoms();
  Delib::ElementTypeCollection elements;
  elements.reserve(N);

  for(AtomIndexType i = 0; i < N; ++i) {
    elements.push_back(_adjacencies[i].elementType);
  }

  return elements;
}

const GraphType& Molecule::getGraph() const {
  return _adjacencies;
}

const StereocenterList& Molecule::getStereocenterList() const {
  return _stereocenters;
}

unsigned Molecule::getNumAdjacencies(const AtomIndexType& a) const {
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

    if(newStereocenter -> numStereopermutations() == 2) {
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
      *this,
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
        stereocenterPtr -> numStereopermutations() == 1
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

  auto activeIndices = getAdjacencies(a);

  temple::inplaceRemoveIf(
    activeIndices,
    [&excludeAdjacent](const auto& adjacentIndex) -> bool {
      return excludeAdjacent.count(adjacentIndex) == 1;
    }
  );

  // Find links between them
  rankingResult.links = GraphAlgorithms::substituentLinks(
    _adjacencies,
    getCycleData(),
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

  auto generateHashes = [](const Molecule& mol) -> std::vector<PseudoHashType> {
    std::vector<PseudoHashType> hashes;

    AtomIndexType N = mol.numAtoms();
    hashes.reserve(N);

    std::vector<BondType> bonds;
    bonds.reserve(Symmetry::constexprProperties::maxSymmetrySize);

    for(AtomIndexType i = 0; i < N; ++i) {
      bonds.clear();

      for(const auto& edge : mol.iterateEdges(i)) {
        bonds.emplace_back(mol.getGraph()[edge].bondType);
      }

      std::sort(
        bonds.begin(),
        bonds.end()
      );

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
        hashAtomEnvironment(
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
  std::unordered_map<PseudoHashType, PseudoHashType> reductionMapping;
  PseudoHashType counter = 0;

  for(const auto& hash : boost::range::join(thisHashes, otherHashes)) {
    if(reductionMapping.count(hash) == 0) {
      reductionMapping.emplace(
        hash,
        counter
      );

      ++counter;
    }
  }

  for(auto& hash : boost::range::join(thisHashes, otherHashes)) {
    hash = reductionMapping.at(hash);
  }

  auto maxHash = counter;

  // This explicit form is needed instead of a lambda for boost's concept checks
  struct HashLookup {
    const std::vector<PseudoHashType>* const hashes;

    using argument_type = AtomIndexType;
    using result_type = PseudoHashType;

    HashLookup(const std::vector<PseudoHashType>& hashes) : hashes(&hashes) {}

    PseudoHashType operator() (const AtomIndexType& i) const {
      return hashes->at(i);
    }
  };

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
        || thisCNSPtr->numStereopermutations() != otherCNSPtr->numStereopermutations()
        || thisCNSPtr->assigned() != otherCNSPtr->assigned()
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
      if(
        stereocenterPtr->numStereopermutations() != otherPtr->numStereopermutations()
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
