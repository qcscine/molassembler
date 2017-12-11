#ifndef INCLUDE_CN_STEREOCENTER_H
#define INCLUDE_CN_STEREOCENTER_H

#include "geometry_assignment/Assignment.h"
#include "symmetry_information/Symmetries.h"
#include "symmetry_information/DynamicProperties.h"

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

namespace Stereocenters {

namespace glue {

RankingInformation::RankedType canonicalize(
  RankingInformation::RankedType sortedSubstituents
);

std::vector<char> makeCanonicalCharacters(
  const RankingInformation::RankedType canonicalizedSubstituents
);

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

boost::optional<const std::vector<unsigned>&> getIndexMapping(
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
  std::vector<AssignmentType> _uniqueAssignmentsCache;

  /*! Mapping between next neighbor atom index to permutational symmetry position
   *
   * Essentially a cached value, valid only for the current assignment
   * (excluding boost::none).
   */
  std::map<AtomIndexType, unsigned> _symmetryPositionMapCache;

public:
/* Public state */
/* Constructors */
  CNStereocenter(
    // The symmetry of this Stereocenter
    const Symmetry::Name& symmetry,
    // The atom this Stereocenter is centered on
    const AtomIndexType& centerAtom,
    // Ranking information of substituents
    const RankingInformation& ranking
  );

/* Modification */
  void addSubstituent(
    const AtomIndexType& newSubstituentIndex,
    const RankingInformation& newRanking,
    const Symmetry::Name& newSymmetry,
    const ChiralStatePreservation& preservationOption
  );

  //! Changes the assignment of the stereocenter
  void assign(const boost::optional<unsigned>& assignment) final;

  /*!
   * The assignment is determined based on three-dimensional positions using
   * angle and chiral distortions from the idealized symmetry.
   */
  void fit(
    const Delib::PositionCollection& positions,
    std::vector<Symmetry::Name> excludeSymmetries = {}
  );

  /*!
   * In case a graph modification changes the ranking of this stereocenter's 
   * substituents, it must be redetermined whether the new configuration is a
   * stereocenter and if so, which one.
   */
  void propagateGraphChange(const RankingInformation& newRanking);

  void propagateVertexRemoval(const AtomIndexType& removedIndex) final;

  void removeSubstituent(
    const AtomIndexType& which,
    const RankingInformation& newRanking,
    const Symmetry::Name& newSymmetry,
    const ChiralStatePreservation& preservationOption
  );

/* Information */
  //! Returns the angle between substituent atoms in the idealized symmetry
  double angle(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const final;

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

  const AtomIndexType& getCentralAtomIndex() const;

  //! Returns the underlying symmetry
  Symmetry::Name getSymmetry() const;

  //! Returns a single-element set containing the central atom
  std::vector<AtomIndexType> involvedAtoms() const final;

  /*!
   * Fetches the number of different assignments possible with the current 
   * substituent ranking and connectivity information. This is also the upper 
   * exclusive bound on the assignment indices that can be used to change the 
   * arrangement of substituents.
   */
  unsigned numAssignments() const final;

  //! If the central symmetry group is changed, we must adapt
  void setSymmetry(const Symmetry::Name& symmetryName);

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
