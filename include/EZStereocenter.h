#ifndef INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H
#define INCLUDE_GRAPH_FEATURES_EZ_STEREOCENTER_H

#include "Stereocenter.h"
#include <boost/optional.hpp>
#include "template_magic/templateMagic.h"

/* TODO 
 */

namespace MoleculeManip {

namespace Stereocenters {

/* Has 2 Assignments: 0 = E, 1 = Z
 */
class EZStereocenter : public Stereocenter {
public:

  /* An EZStereocenter consists of six numbers with the sequence:
   * [CenterAtom, Adjacent of higher rank, Adjacent of lower rank] x 2
   * So e.g. [1, 4, 5, 2, 6, 7], _isEOption=false (Z) is
   *
   *    4       6
   *     \     /
   *      1 = 2
   *     /     \
   *    5       7
   *
   */

private:
/* State */
  std::array<AtomIndexType, 6> _indicesAndRank;
  boost::optional<bool> _isEOption;
  const double _dihedralAngleVariance = 5; // in Degrees

  std::vector<
    std::array<AtomIndexType, 4>
  > _equalPriorityDihedralSequences() const;

  std::vector<
    std::array<AtomIndexType, 4>
  > _differentPriorityDihedralSequences() const;

  std::pair<double, double> _cisLimits() const;

  std::pair<double, double> _transLimits() const;

public:
/* Constructors */
  EZStereocenter() = delete;
  EZStereocenter(
    const AtomIndexType& firstCenter,
    const std::vector<AtomIndexType>& firstCenterRanking,
    const AtomIndexType& secondCenter,
    const std::vector<AtomIndexType>& secondCenterRanking
  ); 

/* Modification */
  virtual void assign(const unsigned& assignment) override final;

/* Information */
  virtual double angle(
    const AtomIndexType& i,
    const AtomIndexType& j,
    const AtomIndexType& k
  ) const override final;

  virtual boost::optional<unsigned> assigned() const override final;

  virtual unsigned numAssignments() const override final;

  virtual std::vector<ChiralityConstraintPrototype> chiralityConstraints() const override final;

  virtual std::vector<DihedralLimits> dihedralLimits() const override final;

  virtual std::string info() const override final;

  virtual std::set<AtomIndexType> involvedAtoms() const override final;

  virtual Type type() const override final;

};

} // eo namespace Stereocenters

} // eo namespace

#endif
