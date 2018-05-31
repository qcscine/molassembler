#ifndef INCLUDE_CN_STEREOCENTER_H
#define INCLUDE_CN_STEREOCENTER_H

#include "stereopermutation/GenerateUniques.h"
#include "chemical_symmetries/DynamicProperties.h"
#include "DistanceGeometry/ValueBounds.h"

#include "Stereocenter.h"
#include "AngstromWrapper.h"

using namespace std::string_literals;

/*! @file
 *
 * Coordinative stereocenter class header file. Permits the storage of
 * particular arrangements of bonded atoms around a central atom and their
 * manipulation.
 *
 * Handles the stereopermutation issue, allowing users to cycle through
 * non-mutually-superimposable arrangements of substituents, here called
 * 'assignments'.
 */

/* Forward declarations */

namespace molassembler {

namespace DistanceGeometry {

class MoleculeSpatialModel;
struct ChiralityConstraint;

} // namespace DistanceGeometry

} // namespace molassembler


namespace molassembler {

namespace Stereocenters {

class CNStereocenter : public Stereocenter {
public:
/* Typedefs */
  using StereopermutationType = stereopermutation::Stereopermutation;

/* Auxiliary classes */
  struct PermutationState {
    //! Stably resorted (by set size) ligands ranking
    RankingInformation::RankedLigandsType canonicalLigands;

    //! Character representation of bonding case
    std::vector<char> symbolicCharacters;

    //! Self-referential representation of links
    stereopermutation::Stereopermutation::LinksSetType selfReferentialLinks;

    //! Mapping from ligand index to modeled ligand plane distance
    std::vector<DistanceGeometry::ValueBounds> ligandDistances;

    using ConeAngleType = std::vector<
      boost::optional<DistanceGeometry::ValueBounds>
    >;
    //! Mapping from ligand index to cone angle optional
    ConeAngleType coneAngles;

    //! Vector of rotationally unique stereopermutations with associated weights
    stereopermutation::StereopermutationsWithWeights permutations;

    //! Vector of whether permutations are feasible or obviously infeasible
    std::vector<unsigned> feasiblePermutations;

    //! Mapping from ligand index to permutational symmetry position
    std::vector<unsigned> symmetryPositionMap;

    PermutationState() = default;
    PermutationState(
      const RankingInformation& ranking,
      const AtomIndexType centerAtom,
      const Symmetry::Name symmetry,
      const GraphType& graph
    );

    /*! Stably re-sort ranked ligand indices in decreasing set size
     *
     * Necessary to avoid treating e.g. AAB and ABB separately, although the
     * resulting assignments are identical.
     *
     * Example: rankedLigands: {5, 8}, {3}, {1, 2, 4}
     * Result: {1, 2, 4}, {5, 8}, {3}
     */
    static RankingInformation::RankedLigandsType canonicalize(
      RankingInformation::RankedLigandsType rankedLigands
    );

    //! Condense ligand ranking information into canonical characters for symbolic computation
    static std::vector<char> transferToSymbolicCharacters(
      const RankingInformation::RankedLigandsType& canonicalLigands
    );

    //! Make ligand-index based ligands self-referential within canonical ligands
    static stereopermutation::Stereopermutation::LinksSetType selfReferentialTransform(
      const RankingInformation::LinksType& rankingLinks,
      const RankingInformation::RankedLigandsType& canonicalLigands
    );

    static std::vector<unsigned> generateLigandToSymmetryPositionMap(
      const stereopermutation::Stereopermutation& assignment,
      const RankingInformation::RankedLigandsType& canonicalLigands
    );

    static std::vector<unsigned> generateSymmetryPositionToLigandMap(
      const stereopermutation::Stereopermutation& assignment,
      const RankingInformation::RankedLigandsType& canonicalLigands
    );

    /*!
     * Generate the character representation of a particular stereopermutation
     * using its map from ?? TODO
     */
    static std::vector<char> makeStereopermutationCharacters(
      const RankingInformation::RankedLigandsType& canonicalLigands,
      const std::vector<char>& canonicalStereopermutationCharacters,
      const std::vector<unsigned>& ligandsAtSymmetryPositions
    );

    /* WARNING: This has to be a copy-initialized optional. Don't change it, unless
     * you want to sift through -fsanitize=address output to find the bug in the
     * optional propagation. It's not safe to return a reference to within a
     * temporary object where this is used.
     */
    static boost::optional<
      std::vector<unsigned>
    > getIndexMapping(
      const Symmetry::properties::SymmetryTransitionGroup& mappingsGroup,
      const ChiralStatePreservation& preservationOption
    );

    static bool isFeasibleStereopermutation(
      const StereopermutationType& assignment,
      const RankingInformation::RankedLigandsType& canonicalLigands,
      const ConeAngleType& coneAngles,
      const RankingInformation& ranking,
      const Symmetry::Name symmetry,
      const GraphType& graph
    );
  };

private:
/* State */
  //! Ranking information of substituents
  RankingInformation _ranking;

  //! Central atom of the Stereocenter
  AtomIndexType _centerAtom;

  //! The symmetry the stereocenter represents
  Symmetry::Name _symmetry;

  //! The current state of assignment (if or not, and if so, which)
  boost::optional<unsigned> _assignmentOption;

  /* Derivative state (cache) */
  PermutationState _cache;

public:
/* Static functions */
/* Constructors */
  CNStereocenter(
    // The base graph
    const GraphType& graph,
    // The symmetry of this Stereocenter
    const Symmetry::Name& symmetry,
    // The atom this Stereocenter is centered on
    const AtomIndexType centerAtom,
    // Ranking information of substituents
    const RankingInformation& ranking
  );

/* Modification */
  /*!
   * Handles the addition of a new substituent to the stereocenter. If the
   * stereocenter contains chiral state, it is attempted to transfer the state
   * into the new assignment space according to the supplied chiral state
   * preservation options
   */
  void addSubstituent(
    const GraphType& graph,
    const AtomIndexType newSubstituentIndex,
    RankingInformation newRanking,
    const Symmetry::Name& newSymmetry,
    const ChiralStatePreservation& preservationOption
  );

  //! Changes the assignment of the stereocenter
  void assign(boost::optional<unsigned> assignment) final;

  //! Assigns the Stereocenter randomly using relative assignment weights
  void assignRandom() final;

  /*!
   * The assignment is determined based on three-dimensional positions using
   * angle and chiral distortions from the idealized symmetry.
   */
  void fit(
    const GraphType& graph,
    const AngstromWrapper& angstromWrapper,
    std::vector<Symmetry::Name> excludeSymmetries = {}
  );

  /*!
   * In case a graph modification changes the ranking of this stereocenter's
   * substituents, it must be redetermined whether the new configuration is a
   * stereocenter and if so, which assignment corresponds to the previous one.
   */
  void propagateGraphChange(
    const GraphType& graph,
    RankingInformation newRanking
  );

  /*!
   * Prepares for the removal of an atom on the graph level, which involves
   * the generation of new atom indices.
   */
  void propagateVertexRemoval(const AtomIndexType removedIndex) final;

  /*!
   * Handles the removal of a substituent from the stereocenter. If the
   * stereocenter carries chiral information, a new assignment can be chosen
   * according to the supplide chiral state preservation option.
   */
  void removeSubstituent(
    const GraphType& graph,
    const AtomIndexType which,
    RankingInformation newRanking,
    const Symmetry::Name& newSymmetry,
    const ChiralStatePreservation& preservationOption
  );

  //! If the central symmetry group is changed, we must adapt
  void setSymmetry(
    const Symmetry::Name symmetryName,
    const GraphType& graph
  );

/* Information */
  //! Returns the angle between substituent ligands in the idealized symmetry
  double angle(
    const unsigned i,
    const unsigned j
  ) const;

  void setModelInformation(
    DistanceGeometry::MoleculeSpatialModel& model,
    const std::function<double(const AtomIndexType)> cycleMultiplierForIndex,
    const double looseningMultiplier
  ) const final;

  /*! Returns the permutation index within the set of possible permutations, if set
   *
   * Returns the (public) information of whether the stereocenter is assigned
   * or not, and if so, which assignment it is.
   */
  boost::optional<unsigned> assigned() const final;

  /*! Returns IOP within the set of symbolic ligand permutations
   *
   * This is different to the assignment. The assignment denotes the index
   * within the set of possible (more specifically, not obviously infeasible)
   * stereopermutations.
   */
  boost::optional<unsigned> indexOfPermutation() const final;

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
  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints() const final;

  //! Returns an information string for diagnostic purposes
  std::string info() const final;

  //! Returns an information string for ranking equality checking purposes
  std::string rankInfo() const;

  //! Returns the underlying ranking
  const RankingInformation& getRanking() const;

  //! Returns the underlying symmetry
  Symmetry::Name getSymmetry() const;

  //! Returns a single-element vector containing the central atom
  std::vector<AtomIndexType> involvedAtoms() const final;

  /*! Returns the number of possible permutations
   *
   * Fetches the number of different assignments possible with the current
   * substituent ranking and connectivity information. This is also the upper
   * exclusive bound on the assignment indices that can be used to change the
   * arrangement of substituents.
   */
  unsigned numAssignments() const final;

  /*! Returns the number of symbolic ligand permutations
   *
   * Fetches the number of permutations determined by symbolic ligand
   * calculation, not considering linking or haptic ligand cones.
   */
  unsigned numStereopermutations() const final;

  /*!
   * Returns Stereocenters::Type::CNStereocenter for resolution of sub-types
   * in homogeneous containers containing pointers to the base class.
   */
  Type type() const final;

/* Operators */
  bool operator == (const CNStereocenter& other) const;
  bool operator < (const CNStereocenter& other) const;
};

} // namespace Stereocenters

} // namespace molassembler

#endif
