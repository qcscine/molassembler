/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Handle arrangements of substituents around an atom-centered symmetry
 *
 * Coordinative stereopermutator class header file. Permits the storage of
 * particular arrangements of bonded atoms around a central atom and their
 * manipulation.
 *
 * Handles the stereopermutation issue, allowing users to cycle through
 * non-mutually-superimposable arrangements of substituents, here called
 * 'assignments'.
 */

#ifndef INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_H
#define INCLUDE_MOLASSEMBLER_ATOM_STEREOPERMUTATOR_H

#include "molassembler/Options.h"
#include "molassembler/OuterGraph.h"

#if __cpp_lib_experimental_propagate_const >= 201505
#define MOLASSEMBLER_ENABLE_PROPAGATE_CONST
#include <experimental/propagate_const>
#endif

using namespace std::string_literals;

namespace Scine {

namespace molassembler {
/* Forward declarations */
struct RankingInformation;
struct PermutationState;

namespace DistanceGeometry {

class SpatialModel;
struct ChiralityConstraint;

} // namespace DistanceGeometry

/**
 * @brief Handles the steric permutation of substituents of a non-terminal
 *   central atom
 *
 * This class handles the permutation of ranked ligands around a central atom.
 * It models its haptic ligands' binding sites and bridges in multidentate
 * ligands in order to decide which stereopermutations are feasible. A
 * stereopermutation may be infeasible, i.e. not realizable in
 * three-dimensional space, if either haptic ligands would intersect due to
 * too close ligand angles for their spatial heft, or if a multidentate
 * ligand's bridge length between binding sites were too short to match the
 * angle. The list of stereopermutations reduced by infeasible
 * stereopermutations is then re-indexed and those indices into the list are
 * called assignments.
 *
 * A Stereopermutator can be unassigned, i.e. the distinct stereopermutation
 * that the substituents are can be indeterminate. If you choose to generate
 * conformers for a molecule that includes unassigned stereopermutators, every
 * conformer will choose an assignment from the pool of feasible assignments
 * randomly, but consistent with relative statistical occurrence weights.
 *
 * E.g. a square planar AABC ligand set will have an A-A cis stereopermutation
 * that occurs twice as often as the A-A trans stereopermutation.
 *
 * @note An instance of this class on a given central atom does not indicate
 *   that that atom is a stereocenter. That is only the case if there are
 *   multiple stereopermutations of the ranked substituents / ligands.
 */
class AtomStereopermutator {
public:
  using PropagatedState = std::tuple<RankingInformation, PermutationState, boost::optional<unsigned>>;

//!@name Special member functions
//!@{
  AtomStereopermutator(AtomStereopermutator&& other) noexcept;
  AtomStereopermutator& operator = (AtomStereopermutator&& other) noexcept;
  AtomStereopermutator(const AtomStereopermutator& other);
  AtomStereopermutator& operator = (const AtomStereopermutator& other);
  ~AtomStereopermutator();

  /*!
   * @brief Construct an AtomStereopermutator
   *
   * @param graph The molecule's graph. This information is needed to model
   *   haptic ligands.
   * @param symmetry The local idealized symmetry to model. Typically the
   *   result of Molecule's determineLocalSymmetry.
   * @param centerAtom The atom index within the molecule that is the center of
   *   the local idealized symmetry
   * @param ranking The ranking of the central atom's substituents and ligand
   *   sites. Typically the result of Molecule's rankPriority.
   */
  AtomStereopermutator(
    // The base graph
    const OuterGraph& graph,
    // The symmetry of this Stereopermutator
    Symmetry::Name symmetry,
    // The atom this Stereopermutator is centered on
    AtomIndex centerAtom,
    // Ranking information of substituents
    RankingInformation ranking
  );
//!@}

//!@name Static functions
//!@{
  /*!
   * @brief Picks a symmetry retaining as much chiral state as possible on a
   *   symmetry position increase
   * @throws std::logic_error If there are no larger symmetries
   */
  static Symmetry::Name up(Symmetry::Name symmetryName);

  /*!
   * @brief Picks a symmetry retaining as much chiral state as possible on a
   *   symmetry position decrease
   * @throws std::logic_error If there are no smaller symmetries
   */
  static Symmetry::Name down(Symmetry::Name symmetryName, unsigned removedSymmetryPosition);
//!@}

//!@name Modifiers
//!@{
  /*!
   * @brief Add a new substituent to the permutator, propagating chiral state
   *
   * Handles the addition of a new substituent to the stereopermutator. If the
   * stereopermutator contains chiral state, it is attempted to transfer the state
   * into the new assignment space according to the supplied chiral state
   * preservation options
   *
   * @param graph The molecule's graph which this permutator helps model.
   * @param newSubstituentIndex The atom index of the new substituent
   * @param newRanking The updated ranking information (after graph addition of
   *   the new substituent)
   * @param newSymmetry The target symmetry of increased size
   * @param preservationOption The behavioral option deciding how chiral state
   *   is propagated.
   */
  void addSubstituent(
    const OuterGraph& graph,
    AtomIndex newSubstituentIndex,
    RankingInformation newRanking,
    boost::optional<Symmetry::Name> newSymmetryOption,
    ChiralStatePreservation preservationOption
  );

  /*!
   * @brief Changes the assignment of the stereopermutator
   *
   * @param assignment The new assignment of the stereopermutator. May be
   *   @p boost::none, which sets the chiral state as indeterminate. Must be
   *   less than the number of assignments if not None.
   */
  void assign(boost::optional<unsigned> assignment);

  /*!
   * @brief Assign the Stereopermutator randomly using relative statistical weights
   *
   * Stereopermutations are generated with relative statistical occurrence
   * weights. The assignment is then chosen from the possible stereopermutations
   * with a discrete distribution whose weights are the corresponding relative
   * statistical occurrences.
   *
   * @note If the stereocenter is already assigned, it is reassigned.
   */
  void assignRandom();

  /**
   * @brief Applies an atom index permutation
   *
   * @param permutation The permutation to apply
   */
  void applyPermutation(const std::vector<AtomIndex>& permutation);

  /*!
   * @brief Determine the symmetry and assignment realized in positions
   *
   * The symmetry and assignment are determined based on three-dimensional
   * positions using angle and chiral distortions from the respective idealized
   * symmetries.
   *
   * @param graph The molecule's graph which this permutator helps model
   * @param angstromWrapper The wrapped positions
   * @param excludeSymmetries Any symmetries that should be excluded from
   *   the fitting procedure
   *
   * @note If Options::tauCriterion is set to @p Enable, this function may
   *   exclude some symmetries from the fitting procedure based on geometric
   *   criteria.
   */
  void fit(
    const OuterGraph& graph,
    const AngstromWrapper& angstromWrapper
  );

  /*!
   * @brief Propagate the stereocenter state through a possible ranking change
   *
   * In case a graph modification changes the ranking of this stereopermutator's
   * substituents, it must be redetermined whether the new configuration is a
   * stereopermutator and if so, which assignment corresponds to the previous one.
   */
  boost::optional<PropagatedState> propagateGraphChange(
    const OuterGraph& graph,
    RankingInformation newRanking
  );

  /*!
   * @brief Adapts atom indices in the internal state to the removal of an atom
   *
   * Atom indices are adapted to a graph-level removal of an atom. The removed
   * index is changed to a placeholder value.
   */
  void propagateVertexRemoval(AtomIndex removedIndex);

  /*!
   * @brief Removes a substituent, propagating state to the new smaller symmetry
   *
   * Handles the removal of a substituent from the stereopermutator. If the
   * stereopermutator carries chiral information, a new assignment can be chosen
   * according to the supplide chiral state preservation option.
   *
   * @warning This should be called after the removal of an atom on the graph level
   */
  void removeSubstituent(
    const OuterGraph& graph,
    AtomIndex which,
    RankingInformation newRanking,
    boost::optional<Symmetry::Name> newSymmetryOption,
    ChiralStatePreservation preservationOption
  );

  /*!
   * @brief Change the symmetry of the permutator
   *
   * @post The permutator is unassigned (chiral state is discarded)
   */
  void setSymmetry(
    Symmetry::Name symmetryName,
    const OuterGraph& graph
  );
//!@}

//!@name Information
//!@{
  /*!
   * @brief Fetches angle between substituent ligands in the idealized symmetry
   *
   * @param i Ligand index one
   * @param j Ligand index two
   *
   * @pre @p i and @p j are valid ligand indices into the underlying
   * RankingInformation's RankingInformation#ligands member.
   *
   * @sa getRanking()
   */
  double angle(unsigned i, unsigned j) const;

  /*!
   * @brief Returns the permutation index within the set of possible permutations, if set
   *
   * Returns the information of whether the stereopermutator is assigned
   * or not, and if so, which assignment it is.
   */
  boost::optional<unsigned> assigned() const;

  //! Returns the central atom this permutator is placed on
  AtomIndex centralIndex() const;

  /*!
   * @brief Returns IOP within the set of symbolic ligand permutations
   *
   * This is different to the assignment. The assignment denotes the index
   * within the set of possible (more specifically, not obviously infeasible)
   * stereopermutations.
   */
  boost::optional<unsigned> indexOfPermutation() const;

  /*!
   * @brief Returns a minimal representation of chirality constraints
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

  /*!
   * @brief Returns the underlying PermutationState
   * @note This is library-internal and not part of the public API
   */
  const PermutationState& getPermutationState() const;

  //! Returns the underlying ranking
  const RankingInformation& getRanking() const;

  //! Returns the underlying symmetry
  Symmetry::Name getSymmetry() const;

  /*!
   * @brief Yields the mapping from ligand indices to symmetry positions
   * @throws std::logic_error if the stereopermutator is unassigned.
   */
  std::vector<unsigned> getSymmetryPositionMap() const;

  /*!
   * @brief Returns the number of possible assignments
   *
   * The number of possible assignments is the number of non-superposable
   * arrangements of the abstract ligand case reduced by trans-arranged
   * multidentate pairs where the bridge length is too short or overlapping
   * haptic cones.
   *
   * For instance, if octahedral M[(A-A)3], there are four abstract arrangements
   * - trans-trans-trans
   * - trans-cis-cis
   * - 2x cis-cis-cis (Δ and Λ isomers, ship propeller-like chirality)
   *
   * However, the number of stereopermutations for a concrete case in which the
   * bridges are too short to allow trans bonding is reduced by all
   * arrangements containing a trans-bonded bidentate ligand, i.e. only Δ and Λ
   * remain. The number of assignments is then only two.
   *
   * This is the upper exclusive bound on Some-type arguments to assign().
   */
  unsigned numAssignments() const;

  /*!
   * @brief Returns the number of possible stereopermutations
   *
   * The number of possible stereopermutations is the number of
   * non-superposable arrangements of the abstract ligand case without removing
   * trans-arranged multidentate pairs or overlapping haptic cones.
   *
   * For instance, if octahedral M[(A-A)3], there are four abstract arrangements
   * - trans-trans-trans
   * - trans-cis-cis
   * - 2x cis-cis-cis (Δ and Λ isomers, ship propeller-like chirality)
   *
   * However, the number of assignments for a concrete case in which the bridges
   * are too short to allow trans binding is reduced by all arrangements
   * containing a trans-bonded bidentate ligand, i.e. only Δ and Λ remain.
   *
   * Fetches the number of permutations determined by symbolic ligand
   * calculation, not considering linking or haptic ligand cones.
   */
  unsigned numStereopermutations() const;
//!@}


//!@name Operators
//!@{
  /**
   * @brief Checks whether the underlying symmetry, central atom index, number
   *   of stereopermutations and current assignment match
   *
   * @param other The other atom stereopermutator to compare against
   *
   * @return Whether the underlying symmetry, central atom index, number of
   *   stereopermutations and current assignment match
   */
  bool operator == (const AtomStereopermutator& other) const;
  //! Inverts operator ==
  bool operator != (const AtomStereopermutator& other) const;

  /**
   * @brief Lexicographically compares the central atom index, the symmetry,
   *   the number of assignments, and the current assignment
   *
   * @param other The other atom stereopermutator to compare against
   *
   * @return Which atom stereopermutator is lexicographically smaller
   */
  bool operator < (const AtomStereopermutator& other) const;
//!@}

private:
  class Impl;

#ifdef MOLASSEMBLER_ENABLE_PROPAGATE_CONST
  std::experimental::propagate_const<
    std::unique_ptr<Impl>
  > _pImpl;
#else
  std::unique_ptr<Impl> _pImpl;
#endif
};

} // namespace molassembler

} // namespace Scine

#endif
