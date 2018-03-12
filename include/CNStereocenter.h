#ifndef INCLUDE_CN_STEREOCENTER_H
#define INCLUDE_CN_STEREOCENTER_H

#include "geometry_assignment/GenerateUniques.h"
#include "symmetry_information/Symmetries.h"
#include "symmetry_information/DynamicProperties.h"
#include "DistanceGeometry/ValueBounds.h"

#include "Stereocenter.h"

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

namespace MoleculeManip {

// Forward-declare Molecule
class Molecule;

namespace Stereocenters {

namespace adhesive {

RankingInformation::RankedType ligandRanking(
  const RankingInformation::RankedType& sortedSubstituents,
  const RankingInformation::RankedType& ligands
);

std::vector<char> canonicalCharacters(
  const RankingInformation::RankedType& rankedLigands
);

UniqueAssignments::Assignment::LinksSetType canonicalLinks(
  const RankingInformation::RankedType& ligands,
  const RankingInformation::RankedType& rankedLigands,
  const RankingInformation::LinksType& rankingLinks
);

} // namespace adhesive

namespace glue {

/*! Stably re-sort ranked substituents in decreasing set size
 *
 * Necessary to avoid treating e.g. AAB and ABB separately, although the
 * resulting assignments are identical.
 *
 * Example: sortedSubstituents: {5, 8}, {3}, {1, 2, 4}
 * Result: {1, 2, 4}, {5, 8}, {3}
 */
RankingInformation::RankedType canonicalize(
  RankingInformation::RankedType sortedSubstituents
);

/*! Transforms canonicalized substituents to character space
 *
 * Example: canonical ranked substituents: {1, 2, 4}, {5, 8}, {3}
 * Result: AAABBC
 */
std::vector<char> makeCanonicalCharacters(
  const RankingInformation::RankedType& canonicalizedSubstituents
);

/*! Transform links from atom-index space to character space
 *
 * In the determination of unique assignments, we have to construct an initial
 * assignment that specifies the ranking differences of ligands and their 
 * links. From the GraphAlgorithm substituentLinks, we receive atom index based
 * link pairs. These need to be altered to character-referential links.
 *
 * Example: canonical ranked substituents: {1, 2, 4}, {5, 8}, {3} == AAABBC
 *          links: {1, 5}, {2, 5}, {4, 5}
 * Result: {0, 3}, {1, 3}, {2, 3} == A1-B1, A2-B1, A3-B1
 */
UniqueAssignments::Assignment::LinksSetType makeLinksSelfReferential(
  const RankingInformation::RankedType& canonicalizedSubstituents,
  const RankingInformation::LinksType& rankingLinks
);

std::map<AtomIndexType, unsigned> makeSymmetryPositionMap(
  const UniqueAssignments::Assignment& assignment,
  const RankingInformation::RankedType& canonicalizedSubstituents
);

std::vector<AtomIndexType> mapToSymmetryPositions(
  const UniqueAssignments::Assignment& assignment,
  const RankingInformation::RankedType& canonicalizedSubstituents
);

std::vector<char> makeAssignmentCharacters(
  const RankingInformation::RankedType& canonicalizedSubstituents,
  const std::vector<char>& canonicalizedAssignmentCharacters,
  const std::vector<AtomIndexType>& atomsAtSymmetryPositions
);

/* WARNING: This has to be a copy-initialized optional. Don't change it, unless
 * you want to sift through -fsanitize=address output to find the bug in the
 * optional propagation. It's not safe to return a reference to within a
 * temporary object where this is used.
 */
boost::optional<std::vector<unsigned>> getIndexMapping(
  const Symmetry::properties::SymmetryTransitionGroup& mappingsGroup,
  const ChiralStatePreservation& preservationOption
);

} // namespace glue

class CNStereocenter : public Stereocenter {
private:
/* Typedefs */
  using AssignmentType = UniqueAssignments::Assignment;

/* Private State */
  //! Ranking information of substituents
  RankingInformation _ranking;

  //! Central atom of the Stereocenter
  AtomIndexType _centerAtom;

  //! The symmetry the stereocenter represents
  Symmetry::Name _symmetry;

  //! The current state of assignment (if or not, and if so, which)
  boost::optional<unsigned> _assignmentOption;

  /* Derivative state (cache) */
  /*! Vector of rotationally unique Assignments
   *
   * Essentially a cached value, valid only for the current ranking.
   */
  UniqueAssignments::AssignmentsWithWeights _uniqueAssignmentsCache;

  /*! Mapping between next neighbor atom index to permutational symmetry position
   *
   * Essentially a cached value, valid only for the current assignment
   * (excluding boost::none).
   */
  std::map<AtomIndexType, unsigned> _symmetryPositionMapCache;

  void _removeImpossibleAssignments(
    UniqueAssignments::AssignmentsWithWeights& data,
    const Molecule& molecule
  );

public:
/* Static functions */
  bool isFeasibleAssignment(
    const AssignmentType& assignment,
    const Molecule& molecule,
    const Symmetry::Name& symmetry,
    const RankingInformation& ranking
  );

/* Constructors */
  CNStereocenter(
    // The base molecule
    const Molecule& molecule,
    // The symmetry of this Stereocenter
    const Symmetry::Name& symmetry,
    // The atom this Stereocenter is centered on
    const AtomIndexType& centerAtom,
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
    const Molecule& molecule,
    const AtomIndexType& newSubstituentIndex,
    const RankingInformation& newRanking,
    const Symmetry::Name& newSymmetry,
    const ChiralStatePreservation& preservationOption
  );

  //! Changes the assignment of the stereocenter
  void assign(const boost::optional<unsigned>& assignment) final;

  //! Assigns the Stereocenter randomly using relative assignment weights
  void assignRandom();

  /*!
   * The assignment is determined based on three-dimensional positions using
   * angle and chiral distortions from the idealized symmetry.
   */
  void fit(
    const Molecule& molecule,
    const Delib::PositionCollection& positions,
    std::vector<Symmetry::Name> excludeSymmetries = {}
  );

  /*!
   * In case a graph modification changes the ranking of this stereocenter's
   * substituents, it must be redetermined whether the new configuration is a
   * stereocenter and if so, which assignment corresponds to the previous one.
   */
  void propagateGraphChange(
    const Molecule& molecule,
    const RankingInformation& newRanking
  );

  /*!
   * Prepares for the removal of an atom on the graph level, which involves
   * the generation of new atom indices.
   */
  void propagateVertexRemoval(const AtomIndexType& removedIndex) final;

  /*!
   * Handles the removal of a substituent from the stereocenter. If the
   * stereocenter carries chiral information, a new assignment can be chosen
   * according to the supplide chiral state preservation option.
   */
  void removeSubstituent(
    const Molecule& molecule,
    const AtomIndexType& which,
    const RankingInformation& newRanking,
    const Symmetry::Name& newSymmetry,
    const ChiralStatePreservation& preservationOption
  );

  //! If the central symmetry group is changed, we must adapt
  void setSymmetry(
    const Molecule& molecule,
    const Symmetry::Name& symmetryName
  );

/* Information */
  //! Returns the angle between substituent atoms in the idealized symmetry
  double angle(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const final;

#if false
  DistanceGeometry::ValueBounds angles(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const;
#endif

  /*!
   * Returns the (public) information of whether the stereocenter is assigned
   * or not, and if so, which assignment it is.
   */
  boost::optional<unsigned> assigned() const final;

  //! Generates a list of chirality constraints on its substituents for DG
  std::vector<ChiralityConstraintPrototype> chiralityConstraints() const final;

  /*!
   * Generates a list of dihedral constraints on its substituents for DG. This
   * is unused here, always returns an empty list.
   */
  std::vector<DihedralLimits> dihedralLimits() const final;

  //! Returns an information string for diagnostic purposes
  std::string info() const final;

  //! Returns an information string for ranking equality checking purposes
  std::string rankInfo() const;

  //! Returns the underlying symmetry
  Symmetry::Name getSymmetry() const;

  //! Returns a single-element vector containing the central atom
  std::vector<AtomIndexType> involvedAtoms() const final;

  /*!
   * Fetches the number of different assignments possible with the current 
   * substituent ranking and connectivity information. This is also the upper 
   * exclusive bound on the assignment indices that can be used to change the 
   * arrangement of substituents.
   */
  unsigned numAssignments() const final;

  /*!
   * Returns Stereocenters::Type::CNStereocenter for resolution of sub-types
   * in homogeneous containers containing pointers to the base class.
   */
  Type type() const final;

/* Operators */
  bool operator == (const CNStereocenter& other) const;
  bool operator < (const CNStereocenter& other) const;

  friend std::basic_ostream<char>& MoleculeManip::operator << (
    std::basic_ostream<char>& os,
    const std::shared_ptr<MoleculeManip::Stereocenters::Stereocenter>& stereocenterPtr
  );
};

} // namespace Stereocenters

} // namespace MoleculeManip

#endif
