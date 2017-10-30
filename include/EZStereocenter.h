#ifndef INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H
#define INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H

#include "Stereocenter.h"

/*! @file
 *
 * Contains the EZStereocenter class declaration, which models E/Z double bond
 * stereocenters in molecules.
 */

/* TODO
 * - Could make dihedralAngleVariance constexpr, but unsure about effects on 
 *   static initialization and if extra declaration required
 */

namespace MoleculeManip {

namespace Stereocenters {

/*! Models stereogenic E/Z conformers in molecular graphs.
 *
 * Has either one or two assignments:
 *
 * - When on either side both substituents have equal rank, there is only one
 *   assignment.
 * - Otherwise, there are two assignments.
 * - In case there is only one substituent on one side, it is considered
 *   "unequal to nothingness". This is so that e.g. diazene::
 *
 *     H
 *       \
 *        N = N
 *              \
 *                H
 *
 *   is properly considered as having E and Z conformers
 *
 * Implementation-wise, _isEOption works as true -> E, false -> Z
 *
 * EZStereocenter is explicitly also instantiable on bonds that are not 
 * stereogenic so that in DG, chirality constraints and dihedral limits that
 * result from the double bond are still collectible.
 */
class EZStereocenter : public Stereocenter {
public:
  using IndexPairsSet = std::set<
    std::pair<AtomIndexType, AtomIndexType>
  >;

  /* An EZStereocenter consists of four to six atom indices:
   *
   * leftCenter,
   * leftHighPriority,
   * leftLowPriority (optional),
   * rightCenter,
   * rightHighPriority,
   * rightLowPriority (optional),
   *
   */

private:
/* State */
  //! Stores the indices of the minimal indices required to have a stereocenter
  AtomIndexType _leftCenter, _leftHighPriority,
                _rightCenter, _rightHighPriority;

  /*! 
   * If there are additional unequal substituents on both sides, their indices
   * graph indices are stored here
   */
  boost::optional<AtomIndexType> _leftLowPriority, _rightLowPriority;

  //! Stores the current assignment state of the stereocenter
  boost::optional<bool> _isEOption;

  //! Stores the number of possible assignments
  unsigned _numAssignments;

  //! Static tolerance in dihedral angle for DG
  static const double _dihedralAngleVariance; // in Degrees

/* Private members */
  //! Generates the dihedral sequences of equal priority (i.e. high with high)
  std::vector<
    std::array<AtomIndexType, 4>
  > _equalPriorityDihedralSequences() const;

  //! Generates the dihedral sequences of unequal priortiy (i.e. high with low)
  std::vector<
    std::array<AtomIndexType, 4>
  > _differentPriorityDihedralSequences() const;

  //! Generates cis dihedral limits
  std::vector<DihedralLimits> _cisDihedralLimits() const;

  //! Generates trans dihedral limits
  std::vector<DihedralLimits> _transDihedralLimits() const;

public:
/* Constructors */
  EZStereocenter() = delete;
  EZStereocenter(
    const AtomIndexType& firstCenter,
    const RankingInformation& firstCenterRanking,
    const AtomIndexType& secondCenter,
    const RankingInformation& secondCenterRanking
  ); 

/* Modification */
  void adaptToRankingChange(const RankingInformation& newRanking) final;

  void assign(const boost::optional<unsigned>& assignment) final;

  void fit(const Delib::PositionCollection& positions) final;

/* Information */
  double angle(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const final;

  boost::optional<unsigned> assigned() const final;

  unsigned numAssignments() const final;

  std::vector<ChiralityConstraintPrototype> chiralityConstraints() const final;

  std::vector<DihedralLimits> dihedralLimits() const final;

  std::string info() const final;

  std::string rankInfo() const;

  std::set<AtomIndexType> involvedAtoms() const final;

  Type type() const final;

/* Operators */
  bool operator == (const EZStereocenter& other) const;
  bool operator < (const EZStereocenter& other) const;

};

} // namespace Stereocenters

} // namespace MoleculeManip

#endif
