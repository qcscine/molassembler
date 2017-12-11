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
 * Implementation-wise, _isZOption works as true -> Z, false -> E
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

private:
/* State */
  //! Stores the indices of the minimal indices required to have a stereocenter
  AtomIndexType _leftCenter, _rightCenter; 

  //! Store the ranking of indices
  RankingInformation _leftRanking, _rightRanking;

  //! Stores the current assignment state of the stereocenter
  boost::optional<bool> _isZOption;

  //! Static tolerance in dihedral angle for DG
  static const double _dihedralAngleVariance; // in Degrees

/* Private members */
  inline AtomIndexType _leftHighPriority() const {
    /* NOTE: high priority is assigned from the back of the ranking since the list
     * is ordered ascending
     *
     * NOTE: back-back always works out to the right selection for all valid
     * inputs.
     */
    return _leftRanking.sortedSubstituents.back().back();
  }

  inline AtomIndexType _leftLowPriority() const {
    return _leftRanking.sortedSubstituents.front().front();
  }

  inline AtomIndexType _rightHighPriority() const {
    return _rightRanking.sortedSubstituents.back().back();
  }

  inline AtomIndexType _rightLowPriority() const {
    return _rightRanking.sortedSubstituents.front().front();
  }

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

  static unsigned _numIndices(const RankingInformation& ranking);

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
  void addSubstituent(
    const AtomIndexType& center,
    const RankingInformation& centerRanking
  );

  void assign(const boost::optional<unsigned>& assignment) final;

  void fit(const Delib::PositionCollection& positions);

  void propagateGraphChange(
    const RankingInformation& firstCenterRanking,
    const RankingInformation& secondCenterRanking
  );

  void propagateVertexRemoval(const AtomIndexType& removedIndex) final;

  void removeSubstituent(
    const AtomIndexType& center,
    const AtomIndexType& which
  );

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

  std::vector<AtomIndexType> involvedAtoms() const final;

  Type type() const final;

/* Operators */
  bool operator == (const EZStereocenter& other) const;
  bool operator < (const EZStereocenter& other) const;

};

} // namespace Stereocenters

} // namespace MoleculeManip

#endif
