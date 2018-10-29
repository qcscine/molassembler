#ifndef INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_IMPL_H
#define INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_IMPL_H

#include "molassembler/BondStereopermutator.h"

#include "boost/optional.hpp"

#include "stereopermutation/Composites.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"

namespace molassembler {

struct BondStereopermutator::Impl {
  Impl(
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    BondIndex edge
  );

  void assign(boost::optional<unsigned> assignment);

  void assignRandom();

  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB
  );

/* Information */
  boost::optional<unsigned> assigned() const;

  const stereopermutation::Composite& composite() const;

  bool hasSameCompositeOrientation(const BondStereopermutator::Impl& other) const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::vector<DistanceGeometry::ChiralityConstraint> chiralityConstraints(
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB
  ) const;

  std::string info() const;

  std::string rankInfo() const;

  BondIndex edge() const;

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
    const AtomStereopermutator& focalStereopermutator,
    const AtomStereopermutator& attachedStereopermutator
  );
};

} // namespace molassembler

#endif
