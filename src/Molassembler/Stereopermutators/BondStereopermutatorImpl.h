/*! @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Private implementation of BondStereopermutator
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_IMPL_H
#define INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_IMPL_H

#include "Molassembler/BondStereopermutator.h"

#include "boost/optional.hpp"

#include "Molassembler/Stereopermutation/Composites.h"
#include "Molassembler/DistanceGeometry/DistanceGeometry.h"

namespace Scine {
namespace Molassembler {

class StereopermutatorList;

struct BondStereopermutator::Impl : public Temple::Crtp::LexicographicComparable<Impl> {
  /* Make sure static_cast-ing our Alignment enum to Composite's Alignment is
   * safe and correct:
   */
  static_assert(
    std::is_same<
      std::underlying_type_t<Stereopermutations::Composite::Alignment>,
      std::underlying_type_t<BondStereopermutator::Alignment>
    >::value,
    "Underlying type of stereopermutation's Alignment and BondStereopermutator's must match"
  );
  static_assert(
    static_cast<Stereopermutations::Composite::Alignment>(BondStereopermutator::Alignment::Staggered) == Stereopermutations::Composite::Alignment::Staggered,
    "Staggered Alignment values do not match across Alignment types"
  );
  static_assert(
    static_cast<Stereopermutations::Composite::Alignment>(BondStereopermutator::Alignment::Eclipsed) == Stereopermutations::Composite::Alignment::Eclipsed,
    "Eclipsed Alignment values do not match across Alignment types"
  );
  static_assert(
    static_cast<Stereopermutations::Composite::Alignment>(BondStereopermutator::Alignment::EclipsedAndStaggered) == Stereopermutations::Composite::Alignment::EclipsedAndStaggered,
    "EclipsedAndStaggered Alignment values do not match across Alignment types"
  );
  static_assert(
    static_cast<Stereopermutations::Composite::Alignment>(BondStereopermutator::Alignment::BetweenEclipsedAndStaggered) == Stereopermutations::Composite::Alignment::BetweenEclipsedAndStaggered,
    "BetweenEclipsedAndStaggered Alignment values do not match across Alignment types"
  );

//!@name Static methods
//!@{
  /*! @brief Check whether a cycle is obviously infeasible
   *
   * @complexity{@math{\Theta(1)}, but not instant, either}
   */
  static bool cycleObviouslyInfeasible(
    const PrivateGraph& graph,
    const StereopermutatorList& stereopermutators,
    const AtomStereopermutator& firstStereopermutator,
    const AtomStereopermutator& secondStereopermutator,
    std::tuple<AtomIndex, AtomIndex, double> dihedral,
    const RankingInformation::Link& link
  );

  static std::vector<unsigned> notObviouslyInfeasibleStereopermutations(
    const PrivateGraph& graph,
    const StereopermutatorList& stereopermutators,
    const Stereopermutations::Composite& composite
  );
//!@}

  /*!
   * @brief Constructor for use of BondStereopermutator in isolation
   *
   * This constructor does not check whether stereopermutations are feasible
   */
  Impl(
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    BondIndex edge,
    Alignment alignment
  );

  /*!
   * @brief Constructor for use of BondStereopermutator in a Molecule
   *
   * This constructor checks whether its stereopermutations are feasible!
   */
  Impl(
    const PrivateGraph& graph,
    const StereopermutatorList& stereopermutators,
    BondIndex edge,
    Alignment alignment
  );

  void assign(boost::optional<unsigned> assignment);

  void assignRandom(Random::Engine& engine);

  void applyPermutation(const std::vector<AtomIndex>& permutation);

  void fit(
    const SitePositionsPair& sitePositions,
    std::pair<FittingReferences, FittingReferences> fittingReferences,
    FittingMode mode
  );

  void propagateGraphChange(
    const AtomStereopermutator::PropagatedState& oldPermutatorState,
    const AtomStereopermutator& newPermutator,
    const PrivateGraph& graph,
    const StereopermutatorList& permutators
  );

  void propagateVertexRemoval(AtomIndex removedIndex);

/* Information */
  Alignment alignment() const;

  boost::optional<unsigned> assigned() const;

  const Stereopermutations::Composite& composite() const;

  double dihedral(
    const AtomStereopermutator& stereopermutatorA,
    SiteIndex siteIndexA,
    const AtomStereopermutator& stereopermutatorB,
    SiteIndex siteIndexB
  ) const;

  bool hasSameCompositeOrientation(const BondStereopermutator::Impl& other) const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::string info() const;

  std::string rankInfo() const;

  BondIndex placement() const;

/* Operators */
  inline auto tie() const {
    return std::make_tuple(std::ref(composite_), assigned());
  }

private:
  /*!
   * Object representing union of both atom stereopermutators, yields abstract
   * spatial arrangements
   */
  Stereopermutations::Composite composite_;
  //! Edge this stereopermutator is placed on
  BondIndex edge_;
  //! List of indices into allPermutations of composite_ that are not obviously infeasible
  std::vector<unsigned> feasiblePermutations_;
  //! Index optional into feasiblePermutations_ of the current assignment
  boost::optional<unsigned> assignment_;

  //! Yields abstract site characters at their shape positions
  static Stereopermutations::Stereopermutation::Occupation makeOccupation_(
    const RankingInformation::RankedSitesType& sitesRanking,
    const AtomStereopermutator::ShapeMap& shapeVertexMap
  );

  static Stereopermutations::Composite constructComposite_(
    const StereopermutatorList& stereopermutators,
    BondIndex edge,
    Alignment alignment
  );

  static Stereopermutations::Composite::OrientationState makeOrientationState_(
    const AtomStereopermutator& focalStereopermutator,
    const AtomStereopermutator::ShapeMap& focalShapeMap,
    const AtomStereopermutator& attachedStereopermutator
  );
};

} // namespace Molassembler
} // namespace Scine

#endif
