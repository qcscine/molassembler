#ifndef INCLUDE_MOLASSEMBLER_ATOM_STEREOCENTER_IMPL_H
#define INCLUDE_MOLASSEMBLER_ATOM_STEREOCENTER_IMPL_H

#include "molassembler/AtomStereocenter.h"

#include "molassembler/DistanceGeometry/DistanceGeometry.h"
#include "molassembler/Stereocenters/PermutationState.h"
#include "molassembler/RankingInformation.h"
#include "chemical_symmetries/Names.h"
#include "boost/optional.hpp"

namespace molassembler {

/* forward-declarations */

class AtomStereocenter::Impl {
public:
/* Typedefs */
  using StereopermutationType = stereopermutation::Stereopermutation;

/* Static functions */
/* Constructors */
  Impl(
    // The base graph
    const OuterGraph& graph,
    // The symmetry of this Stereocenter
    Symmetry::Name symmetry,
    // The atom this Stereocenter is centered on
    AtomIndex centerAtom,
    // Ranking information of substituents
    RankingInformation ranking
  );

/* Modification */
  /*!
   * Handles the addition of a new substituent to the stereocenter. If the
   * stereocenter contains chiral state, it is attempted to transfer the state
   * into the new assignment space according to the supplied chiral state
   * preservation options
   */
  void addSubstituent(
    const OuterGraph& graph,
    AtomIndex newSubstituentIndex,
    RankingInformation newRanking,
    Symmetry::Name newSymmetry,
    ChiralStatePreservation preservationOption
  );

  //! Changes the assignment of the stereocenter
  void assign(boost::optional<unsigned> assignment);

  //! Assigns the Stereocenter randomly using relative assignment weights
  void assignRandom();

  /*!
   * The symmetry and assignment are determined based on three-dimensional
   * positions using angle and chiral distortions from the respective idealized
   * symmetries.
   */
  void fit(
    const OuterGraph& graph,
    const AngstromWrapper& angstromWrapper,
    const std::vector<Symmetry::Name>& excludeSymmetries = {}
  );

  /*!
   * In case a graph modification changes the ranking of this stereocenter's
   * substituents, it must be redetermined whether the new configuration is a
   * stereocenter and if so, which assignment corresponds to the previous one.
   */
  void propagateGraphChange(
    const OuterGraph& graph,
    RankingInformation newRanking
  );

  /*!
   * Prepares for the removal of an atom on the graph level, which involves
   * the generation of new atom indices.
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  /*!
   * Handles the removal of a substituent from the stereocenter. If the
   * stereocenter carries chiral information, a new assignment can be chosen
   * according to the supplide chiral state preservation option.
   */
  void removeSubstituent(
    const OuterGraph& graph,
    AtomIndex which,
    RankingInformation newRanking,
    Symmetry::Name newSymmetry,
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
   * Returns the (public) information of whether the stereocenter is assigned
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

  //! Generates a list of chirality constraints on its substituents for DG
  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints(
    double looseningMultiplier
  ) const;

  //! Returns an information string for diagnostic purposes
  std::string info() const;

  //! Returns an information string for ranking equality checking purposes
  std::string rankInfo() const;

  //! Returns the underlying ranking
  const RankingInformation& getRanking() const;

  //! Returns the underlying symmetry
  Symmetry::Name getSymmetry() const;

  /*! Yields the mapping from ligand indices to symmetry positions
   *
   * \throws std::logic_error if the stereocenter is unassigned.
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

  void setModelInformation(
    DistanceGeometry::SpatialModel& model,
    const std::function<double(const AtomIndex)>& cycleMultiplierForIndex,
    double looseningMultiplier
  ) const;


/* Operators */
  bool operator == (const Impl& other) const;
  bool operator < (const Impl& other) const;

private:
/* State */
  //! Ranking information of substituents
  RankingInformation _ranking;

  //! Central atom of the Stereocenter
  AtomIndex _centerAtom;

  //! The symmetry the stereocenter represents
  Symmetry::Name _symmetry;

  //! The current state of assignment (if or not, and if so, which)
  boost::optional<unsigned> _assignmentOption;

  /* Derivative state (cache) */
  PermutationState _cache;
};

} // namespace molassembler

#endif
