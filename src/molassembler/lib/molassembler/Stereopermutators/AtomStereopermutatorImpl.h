/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Private implementation of AtomStereopermutator
 */

#ifndef INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_IMPL_H
#define INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_IMPL_H

#include "molassembler/AtomStereopermutator.h"

#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/Stereopermutators/PermutationState.h"

namespace Scine {

namespace molassembler {

/* forward-declarations */

class AtomStereopermutator::Impl {
public:
/* Typedefs */
  using StereopermutationType = stereopermutation::Stereopermutation;

//!@name Static functions
//!@{
  /*!
   * @brief Picks a symmetry retaining as much chiral state as possible on a
   *   symmetry position increase
   * @throws std::logic_error If there are no larger symmetries
   * @note Behavior is dependent on ChiralStatePreservation option
   */
  static Symmetry::Name up(Symmetry::Name symmetryName);

  /*!
   * @brief Picks a symmetry retaining as much chiral state as possible on a
   *   symmetry position decrease
   * @throws std::logic_error If there are no smaller symmetries
   * @note Behavior is dependent on ChiralStatePreservation option
   */
  static Symmetry::Name down(Symmetry::Name symmetryName, unsigned removedSymmetryPosition);
//!@}

/* Constructors */
  Impl(
    // The base graph
    const OuterGraph& graph,
    // The symmetry of this Stereopermutator
    Symmetry::Name symmetry,
    // The atom this Stereopermutator is centered on
    AtomIndex centerAtom,
    // Ranking information of substituents
    RankingInformation ranking
  );

/* Modification */
  /*!
   * Handles the addition of a new substituent to the stereopermutator. If the
   * stereopermutator contains chiral state, it is attempted to transfer the state
   * into the new assignment space according to the supplied chiral state
   * preservation options
   */
  void addSubstituent(
    const OuterGraph& graph,
    AtomIndex newSubstituentIndex,
    RankingInformation newRanking,
    boost::optional<Symmetry::Name> newSymmetryOption,
    ChiralStatePreservation preservationOption
  );

  //! Changes the assignment of the stereopermutator
  void assign(boost::optional<unsigned> assignment);

  //! Assigns the Stereopermutator randomly using relative assignment weights
  void assignRandom();

  /**
   * @brief Applies an atom index permutation
   *
   * @param permutation The permutation to apply
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  /*!
   * The symmetry and assignment are determined based on three-dimensional
   * positions using angle and chiral distortions from the respective idealized
   * symmetries.
   */
  void fit(
    const OuterGraph& graph,
    const AngstromWrapper& angstromWrapper
  );

  /*!
   * In case a graph modification changes the ranking of this stereopermutator's
   * substituents, it must be redetermined whether the new configuration is a
   * stereopermutator and if so, which assignment corresponds to the previous one.
   */
  boost::optional<PropagatedState> propagateGraphChange(
    const OuterGraph& graph,
    RankingInformation newRanking
  );

  /*!
   * Prepares for the removal of an atom on the graph level, which involves
   * the generation of new atom indices.
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  /*!
   * Handles the removal of a substituent from the stereopermutator. If the
   * stereopermutator carries chiral information, a new assignment can be chosen
   * according to the supplide chiral state preservation option.
   */
  void removeSubstituent(
    const OuterGraph& graph,
    AtomIndex which,
    RankingInformation newRanking,
    boost::optional<Symmetry::Name> newSymmetryOption,
    ChiralStatePreservation preservationOption
  );

  //! If the central symmetry group is changed, we must adapt
  void setSymmetry(
    Symmetry::Name symmetryName,
    const OuterGraph& graph
  );

/* Information */
  //! Returns the angle between substituent ligands in the idealized symmetry
  double angle(unsigned i, unsigned j) const;

  /*! Returns the permutation index within the set of possible permutations, if set
   *
   * Returns the (public) information of whether the stereopermutator is assigned
   * or not, and if so, which assignment it is.
   */
  boost::optional<unsigned> assigned() const;

  //! Returns a single-element vector containing the central atom
  AtomIndex centralIndex() const;

  /*! Returns IOP within the set of symbolic ligand permutations
   *
   * This is different to the assignment. The assignment denotes the index
   * within the set of possible (more specifically, not obviously infeasible)
   * stereopermutations.
   */
  boost::optional<unsigned> indexOfPermutation() const;

  /*! Returns a minimal representation of chirality constraints
   *
   * Every minimal representation consists only of ligand indices.
   *
   * The minimal representation assumes that all Symmetry tetrahedron
   * definitions are defined to be Positive targets, which is checked in
   * the chemical_symmetries tests.
   */
  std::vector<
    std::array<boost::optional<unsigned>, 4>
  > minimalChiralityConstraints() const;

  //! Returns an information string for diagnostic purposes
  std::string info() const;

  //! Returns an information string for ranking equality checking purposes
  std::string rankInfo() const;

  /*! Returns the underlying PermutationState
   *
   * @note This is library-internal and not part of the public API
   */
  const PermutationState& getPermutationState() const;

  //! Returns the underlying ranking
  const RankingInformation& getRanking() const;

  //! Returns the underlying symmetry
  Symmetry::Name getSymmetry() const;

  /*! Yields the mapping from ligand indices to symmetry positions
   *
   * \throws std::logic_error if the stereopermutator is unassigned.
   */
  std::vector<unsigned> getSymmetryPositionMap() const;

  /*! Returns the number of possible permutations
   *
   * Fetches the number of different assignments possible with the current
   * substituent ranking and connectivity information. This is also the upper
   * exclusive bound on the assignment indices that can be used to change the
   * arrangement of substituents.
   */
  unsigned numAssignments() const;

  /*! Returns the number of symbolic ligand permutations
   *
   * Fetches the number of permutations determined by symbolic ligand
   * calculation, not considering linking or haptic ligand cones.
   */
  unsigned numStereopermutations() const;


/* Operators */
  bool operator == (const Impl& other) const;
  bool operator < (const Impl& other) const;

private:
/* State */
  //! Ranking information of substituents
  RankingInformation _ranking;

  //! Central atom of the Stereopermutator
  AtomIndex _centerAtom;

  //! The symmetry the stereopermutator represents
  Symmetry::Name _symmetry;

  //! The current state of assignment (if or not, and if so, which)
  boost::optional<unsigned> _assignmentOption;

  /* Derivative state (cache) */
  PermutationState _cache;
};

} // namespace molassembler

} // namespace Scine

#endif
