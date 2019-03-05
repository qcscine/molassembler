/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Directed conformer generation class and helper functions
 */

#ifndef INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H
#define INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H

#include "molassembler/Types.h"
#include "Utils/Typenames.h"
#include "boost/optional/optional_fwd.hpp"
#include "boost_outcome/outcome.hpp"
#include "boost/variant/variant_fwd.hpp"
#include <map>
#include <memory>
#include <vector>

namespace Scine {
namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

class Molecule;
class BondStereopermutator;

/**
 * @brief
 *
 * @note This type is not copyable.
 */
class DirectedConformerGenerator {
public:
//!@name Public types
//!@{
  using BondList = std::vector<BondIndex>;
  using DecisionList = std::vector<std::uint8_t>;

  /**
   * @brief Reason why a bond is ignored
   */
  enum class IgnoreReason {
    //! There is not an assigned stereopermutator on both ends of the bond
    AtomStereopermutatorPreconditionsUnmet,
    //! There is already an assigned bond stereopermutator on the bond
    HasAssignedBondStereopermutator,
    //! At least one constituting atom is terminal
    HasTerminalConstitutingAtom,
    //! This bond is in a cycle of size three or four
    InSmallCycle,
    //! This bond is an eta bond (indicates bonding to haptic ligands, and therefore excluded)
    IsEtaBond,
    //! Rotation about this bond is isotropic (all ligands have same ranking on at least one side)
    RotationIsIsotropic
  };
//!@}

//!@name Static functions
//!@{
  /**
   * @brief Decide whether to consider a bond's dihedral values for directed
   *   conformer generation or not
   *
   * @param bondIndex The bond to consider
   * @param molecule The molecule in which @p bond exists
   * @param smallestCycleMap A map of atom indices to the smallest cycle
   *   they are in
   *
   * @see makeSmallestCycleMap
   *
   * @return Either a reason why the bond was ignored, or a BondStereopermutator
   *   placed on the suggested bond indicating that the bond should be
   *   considered.
   */
  static boost::variant<IgnoreReason, BondStereopermutator> considerBond(
    const BondIndex& bondIndex,
    const Molecule& molecule,
    const std::map<AtomIndex, unsigned>& smallestCycleMap
  );

  /**
   * @brief Calculates a distance metric between two decision lists for
   *   dihedral permutations
   *
   * @param a The first distance metric
   * @param b The second distance metric
   * @param bounds Upper exclusive bound on values at each position
   *
   * @return A distance metric between @p a and @p b.
   */
  static unsigned distance(
    const DecisionList& a,
    const DecisionList& b,
    const DecisionList& bounds
  );
//!@}

//!@name Constructors
//!@{
  explicit DirectedConformerGenerator(Molecule molecule);
//!@}

//!@name Special member functions
//!@{
  DirectedConformerGenerator(DirectedConformerGenerator&& other) noexcept;
  DirectedConformerGenerator& operator = (DirectedConformerGenerator&& other) noexcept;
  DirectedConformerGenerator(const DirectedConformerGenerator& other) = delete;
  DirectedConformerGenerator& operator = (const DirectedConformerGenerator& other) = delete;
  ~DirectedConformerGenerator();
//!@}

//!@name Modification
//!@{
  /*!
   * @brief Generate a new list of discrete dihedral arrangement choices
   * @throws std::logic_error If the underlying set is full, i.e. all decision
   *   lists for conformers have been generated.
   * @post The new DecisionList is part of the stored list of generated
   *   decision lists and will not be generated again.
   */
  DecisionList generateNewDecisionList();

  /*!
   * @brief Adds a decision list to the underlying set-like data structure
   * @returns @p true if @p decisionList wasn't already part of the set
   */
  bool insert(const DecisionList& decisionList);

  //! Checks whether a DecisionList is part of the underlying set
  bool contains(const DecisionList& decisionList);
//!@}

//!@name Information
//!@{
  //! Accessor for list of relevant bonds, O(1)
  const BondList& bondList() const;

  //! Number of conformer decision lists stored in the underlying set-like data structure
  unsigned conformerCount() const;

  //! Number of conformers needed for full ensemble, O(1)
  unsigned idealEnsembleSize() const;

  //! Try to generate a conformer for a particular decision list
  outcome::result<Utils::PositionCollection> generateConformer(const DecisionList& decisionList);

  //! Infer a decision list from positional information
  DecisionList getDecisionList(Utils::PositionCollection positions) const;
//!@}

private:
  class Impl;
  std::unique_ptr<Impl> _pImpl;
};

} // namespace molassembler
} // namespace Scine

#endif
