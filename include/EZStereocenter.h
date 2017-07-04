#ifndef INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H
#define INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H

#include "Stereocenter.h"
#include <boost/optional.hpp>
#include "template_magic/TemplateMagic.h"

namespace MoleculeManip {

namespace Stereocenters {

/* Has either one or two assignments:
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
 *   can be properly considered as having E and Z conformers
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
  AtomIndexType _leftCenter, _leftHighPriority,
                _rightCenter, _rightHighPriority;
  boost::optional<AtomIndexType> _leftLowPriority, _rightLowPriority;
  boost::optional<bool> _isEOption;
  unsigned _numAssignments;

  static const double _dihedralAngleVariance; // in Degrees

/* Private members */
  std::vector<
    std::array<AtomIndexType, 4>
  > _equalPriorityDihedralSequences() const;

  std::vector<
    std::array<AtomIndexType, 4>
  > _differentPriorityDihedralSequences() const;

  std::vector<DihedralLimits> _cisDihedralLimits() const;
  std::vector<DihedralLimits> _transDihedralLimits() const;

public:
/* Constructors */
  EZStereocenter() = delete;
  EZStereocenter(
    const AtomIndexType& firstCenter,
    const std::vector<AtomIndexType>& firstCenterRanking,
    const IndexPairsSet& firstCenterEqualPairs,
    const AtomIndexType& secondCenter,
    const std::vector<AtomIndexType>& secondCenterRanking,
    const IndexPairsSet& secondCenterEqualPairs
  ); 

/* Modification */
  void assign(const unsigned& assignment) final;

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

  std::set<AtomIndexType> involvedAtoms() const final;

  Type type() const final;

/* Operators */
  bool operator == (const EZStereocenter& other) const;

};

} // namespace Stereocenters

} // namespace MoleculeManip

#endif
