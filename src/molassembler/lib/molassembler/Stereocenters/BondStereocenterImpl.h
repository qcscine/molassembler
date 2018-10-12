#ifndef INCLUDE_MOLASSEMBLER_BOND_STEREOCENTER_IMPL_H
#define INCLUDE_MOLASSEMBLER_BOND_STEREOCENTER_IMPL_H

#include "molassembler/BondStereocenter.h"

#include "boost/optional.hpp"

#include "stereopermutation/Composites.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"

namespace molassembler {

struct BondStereocenter::Impl {
  Impl(
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB,
    BondIndex edge
  );

  void assign(boost::optional<unsigned> assignment);

  void assignRandom();

  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB
  );

/* Information */
  boost::optional<unsigned> assigned() const;

  bool hasSameCompositeOrientation(const BondStereocenter::Impl& other) const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints(
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB
  ) const;

  std::string info() const;

  std::string rankInfo() const;

  BondIndex edge() const;

  void setModelInformation(
    DistanceGeometry::SpatialModel& model,
    const AtomStereocenter& stereocenterA,
    const AtomStereocenter& stereocenterB,
    double looseningMultiplier
  ) const;

/* Operators */
  bool operator == (const Impl& other) const;
  bool operator != (const Impl& other) const;

private:
  stereopermutation::Composite _composite;
  BondIndex _edge;
  boost::optional<unsigned> _assignment;

  //! Yields abstract ligand characters at their symmetry positions
  static std::vector<char> _charifyRankedLigands(
    const std::vector<std::vector<unsigned>>& ligandsRanking,
    const std::vector<unsigned>& symmetryPositionMap
  );

  static stereopermutation::Composite::OrientationState _makeOrientationState(
    const AtomStereocenter& focalStereocenter,
    const AtomStereocenter& attachedStereocenter
  );
};

} // namespace molassembler

#endif
