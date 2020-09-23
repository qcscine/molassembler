/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Directed conformer generation class and helper functions implementation file
 */

#ifndef INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_IMPL_H
#define INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_IMPL_H

#include "Molassembler/DirectedConformerGenerator.h"
#include "Molassembler/Molecule.h"

#include "Molassembler/Temple/BoundedNodeTrie.h"

namespace Scine {
namespace Molassembler {

class DirectedConformerGenerator::Impl {
public:
  static unsigned distance(
    const DecisionList& a,
    const DecisionList& b,
    const DecisionList& bounds
  );

  static boost::variant<IgnoreReason, BondStereopermutator> considerBond(
    const BondIndex& bondIndex,
    const Molecule& molecule,
    const std::unordered_map<AtomIndex, unsigned>& smallestCycleMap,
    BondStereopermutator::Alignment alignment
  );

  Impl(
    Molecule molecule,
    const BondStereopermutator::Alignment alignment,
    const BondList& bondsToConsider
  );

  DecisionList generateNewDecisionList(Random::Engine& engine);

  bool insert(const DecisionList& decisionList) {
    return decisionLists_.insert(decisionList);
  }

  void clear() {
    return decisionLists_.clear();
  }

  bool contains(const DecisionList& decisionList) const {
    return decisionLists_.contains(decisionList);
  }

  const BondList& bondList() const {
    return relevantBonds_;
  }

  unsigned decisionListSetSize() const {
    return decisionLists_.size();
  }

  unsigned idealEnsembleSize() const {
    return decisionLists_.capacity();
  }

  outcome::result<Utils::PositionCollection> checkGeneratedConformation(
    outcome::result<Utils::PositionCollection> conformerResult,
    const DecisionList& decisionList,
    BondStereopermutator::FittingMode fitting
  ) const;

  outcome::result<Utils::PositionCollection> generateRandomConformation(
    const DecisionList& decisionList,
    const DistanceGeometry::Configuration& configuration,
    BondStereopermutator::FittingMode fitting
  ) const;

  outcome::result<Utils::PositionCollection> generateConformation(
    const DecisionList& decisionList,
    unsigned seed,
    const DistanceGeometry::Configuration& configuration,
    BondStereopermutator::FittingMode fitting
  ) const;

  DecisionList getDecisionList(
    const Utils::AtomCollection& atomCollection,
    BondStereopermutator::FittingMode fitting
  ) const;

  DecisionList getDecisionList(
    const Utils::PositionCollection& positions,
    BondStereopermutator::FittingMode fitting
  ) const;

  Molecule conformationMolecule(const DecisionList& decisionList) const;

  void enumerate(
    std::function<void(const DecisionList&, Utils::PositionCollection)> callback,
    unsigned seed,
    const EnumerationSettings& settings
  );

  Relabeler relabeler() const;

  std::vector<int> binMidpointIntegers(const DecisionList& decision) const;

private:
  Molecule molecule_;
  BondStereopermutator::Alignment alignment_;
  BondList relevantBonds_;

  /* This data structure is primarily designed for use as a set-like type that
   * can contain all choices of discrete enumerated dihedral positions for a
   * dihedral chain.
   *
   * Say you have three different dihedrals, each of which may have a distinct
   * number of possible rotations (although in organic molecules, the most
   * common one will naturally be three). Then the bounds for construction of
   * this tree are e.g. {3, 2, 4}.
   *
   * Possible DecisionLists are then e.g.:
   * - {0, 0, 0}, (this is the minimal choice list if ordered lexicographically)
   * - {2, 1, 3}, (this is the maximal choice list if ordered lexicographically)
   * - {1, 1, 3},
   * - ...
   *
   * This data structure can help you keep track of which choices at each
   * dihedral you have explored and which ones might lead to a conformer that
   * is most different from the ones you already have.
   */
  Temple::BoundedNodeTrie<std::uint8_t> decisionLists_;
};

} // namespace Molassembler
} // namespace Scine

#endif
