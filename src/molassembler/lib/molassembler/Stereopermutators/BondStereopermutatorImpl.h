/*! @file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Private implementation of BondStereopermutator
 */

#ifndef INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_IMPL_H
#define INCLUDE_MOLASSEMBLER_BOND_STEREOPERMUTATOR_IMPL_H

#include "molassembler/BondStereopermutator.h"

#include "boost/optional.hpp"

#include "temple/OperatorSuppliers.h"

#include "stereopermutation/Composites.h"
#include "molassembler/DistanceGeometry/DistanceGeometry.h"

namespace Scine {

namespace molassembler {

class StereopermutatorList;
struct LinkInformation;

struct BondStereopermutator::Impl : public temple::crtp::LexicographicComparable<Impl> {
  /* Make sure static_cast-ing our Alignment enum to Composite's Alignment is
   * safe and correct:
   */
  static_assert(
    std::is_same<
      std::underlying_type_t<stereopermutation::Composite::Alignment>,
      std::underlying_type_t<BondStereopermutator::Alignment>
    >::value,
    "Underlying type of stereopermutation's Alignment and BondStereopermutator's must match"
  );
  static_assert(
    static_cast<stereopermutation::Composite::Alignment>(BondStereopermutator::Alignment::Staggered) == stereopermutation::Composite::Alignment::Staggered,
    "Staggered Alignment values do not match across Alignment types"
  );
  static_assert(
    static_cast<stereopermutation::Composite::Alignment>(BondStereopermutator::Alignment::Eclipsed) == stereopermutation::Composite::Alignment::Eclipsed,
    "Eclipsed Alignment values do not match across Alignment types"
  );

//!@name Static methods
//!@{
  /*! @brief Check whether a cycle is obviously infeasible
   *
   * @complexity{@math{\Theta(1)}, but not instant, either}
   */
  static bool cycleObviouslyInfeasible(
    const InnerGraph& graph,
    const StereopermutatorList& stereopermutators,
    const AtomStereopermutator& firstStereopermutator,
    const AtomStereopermutator& secondStereopermutator,
    std::tuple<AtomIndex, AtomIndex, double> dihedral,
    const LinkInformation& link
  );

  static std::vector<unsigned> notObviouslyInfeasibleStereopermutations(
    const InnerGraph& graph,
    const StereopermutatorList& stereopermutators,
    const stereopermutation::Composite& composite
  );
//!@}

  /*!
   * @brief Constructor for use of BondStereopermutator in isolation
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
    const InnerGraph& graph,
    const StereopermutatorList& stereopermutators,
    BondIndex edge,
    Alignment alignment
  );

  void assign(boost::optional<unsigned> assignment);

  void assignRandom(random::Engine& engine);

  void applyPermutation(const std::vector<AtomIndex>& permutation);

  void fit(
    const AngstromWrapper& angstromWrapper,
    const AtomStereopermutator& stereopermutatorA,
    const AtomStereopermutator& stereopermutatorB,
    FittingMode mode
  );

  void propagateGraphChange(
    const AtomStereopermutatorPropagatedState& oldPermutatorState,
    const AtomStereopermutator& newPermutator,
    const InnerGraph& graph,
    const StereopermutatorList& permutators
  );

/* Information */
  Alignment alignment() const;

  boost::optional<unsigned> assigned() const;

  const stereopermutation::Composite& composite() const;

  double dihedral(
    const AtomStereopermutator& stereopermutatorA,
    unsigned siteIndexA,
    const AtomStereopermutator& stereopermutatorB,
    unsigned siteIndexB
  ) const;

  bool hasSameCompositeOrientation(const BondStereopermutator::Impl& other) const;

  boost::optional<unsigned> indexOfPermutation() const;

  unsigned numAssignments() const;

  unsigned numStereopermutations() const;

  std::string info() const;

  std::string rankInfo() const;

  BondIndex edge() const;

/* Operators */
  inline auto tie() const {
    return std::make_tuple(std::ref(_composite), assigned());
  }

private:
  /*!
   * Object representing union of both atom stereopermutators, yields abstract
   * spatial arrangements
   */
  stereopermutation::Composite _composite;
  //! Edge this stereopermutator is placed on
  BondIndex _edge;
  //! List of stereopermutations of _composite that are not obviously infeasible
  std::vector<unsigned> _feasiblePermutations;
  //! Index optional into _feasiblePermutations of the current assignment
  boost::optional<unsigned> _assignment;

  //! Yields abstract site characters at their symmetry positions
  static std::vector<char> _charifyRankedSites(
    const std::vector<std::vector<unsigned>>& sitesRanking,
    const std::vector<unsigned>& symmetryPositionMap
  );

  static stereopermutation::Composite _constructComposite(
    const StereopermutatorList& stereopermutators,
    BondIndex edge,
    Alignment alignment
  );

  static stereopermutation::Composite::OrientationState _makeOrientationState(
    const AtomStereopermutator& focalStereopermutator,
    const AtomStereopermutator& attachedStereopermutator
  );
};

} // namespace molassembler

} // namespace Scine

#endif
