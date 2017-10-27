#ifndef INCLUDE_CN_STEREOCENTER_H
#define INCLUDE_CN_STEREOCENTER_H

#include "geometry_assignment/Assignment.h"
#include "symmetry_information/Symmetries.h"

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

/* TODO
 * - Consider altering _neighborSymmetryPositionMap and _neighborCharMap into
 *   plain vectors of length Symmetry::size. Which way do lookups go and how
 *   would speed be affected?
 * - Some of the private members could be static
 */

class CNStereocenter : public Stereocenter {
private:
/* Typedefs */
  using AssignmentType = UniqueAssignments::Assignment;

/* Private State */
  //! Mapping between next neighbor atom index and symbolic ligand character
  std::map<AtomIndexType, char> _neighborCharMap;
  //! Links between sorted indices
  UniqueAssignments::Assignment::LinksSetType _links;

  //! Mapping between next neighbor atom index to permutational symmetry position
  std::map<AtomIndexType, unsigned> _neighborSymmetryPositionMap;
  //! Vector of rotationally unique Assignments
  std::vector<AssignmentType> _uniqueAssignments;
  
/* Private members */
  /*!
   * Generates the mapping of next neighbor atom graph indices to symmetry
   * positions for a particular assignment
   */
  static std::map<AtomIndexType, unsigned> _makeNeighborSymmetryPositionMap(
    const UniqueAssignments::Assignment& assignment,
    const std::map<AtomIndexType, char> neighborCharMap
  );

  /*!
   * Reduce substituent atoms at central atom to a mapping of their indices to
   * a symbolic ligand character. e.g.
   * map{4 -> 'C', 6 -> 'A', 9 -> 'A', 14 -> 'B'}
   */
  static std::map<AtomIndexType, char> _reduceSubstituents(
    const RankingInformation& ranking
  );

  /*!
   * Reduce a mapping of atom indices to symbolic ligand characters to a 
   * character vector
   */
  std::vector<char> _reduceNeighborCharMap(
    const std::map<AtomIndexType, char>& neighborCharMap
  ) const;

  UniqueAssignments::Assignment::LinksSetType _makeLinks(
    const RankingInformation& ranking
    /*const std::vector<AtomIndexType>& sortedIndices,
    const std::set<
      std::pair<AtomIndexType, AtomIndexType>
    >& linkedPairsSet*/
  ) const;

public:
/* Public state */
  //! The symmetry of the overall stereocenter
  Symmetry::Name symmetry;
  //! Central atom of the Stereocenter, const on assignment
  AtomIndexType centerAtom; 
  //! The current state of assignment (if or not, and if so, which)
  boost::optional<unsigned> assignmentOption;

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
  //! Changes the assignment of the stereocenter
  void assign(const unsigned& assignment) final;

  /*!
   * In case a graph modification changes the ranking of this stereocenter's 
   * substituents, it must be redetermined whether the new configuration is a
   * stereocenter and if so, which one.
   */
  void adaptToRankingChange(const RankingInformation& newRanking) final;

  //! If the central symmetry group is changed, we must adapt
  void changeSymmetry(const Symmetry::Name& symmetryName);

  /*!
   * The assignment is determined based on three-dimensional positions using
   * angle and chiral distortions from the idealized symmetry.
   */
  void fit(const Delib::PositionCollection& positions) final;

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

  //! Returns a single-element set containing the central atom
  std::set<AtomIndexType> involvedAtoms() const final;

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
