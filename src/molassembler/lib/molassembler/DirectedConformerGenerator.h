/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Directed conformer generation class and helper functions
 */

#ifndef INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H
#define INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H

#include "molassembler/Types.h"
#include "Utils/Typenames.h"
#include "boost_outcome/outcome.hpp"
#include <map>
#include <memory>
#include <vector>

namespace Scine {
namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

class Molecule;

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
   * @brief Importance of a bond
   */
  enum class BondImportance {
    Keep,
    IgnoreHasAssignedBondStereopermutator,
    IgnoreInSmallCycle,
    IgnoreRotationIsIsotropic
  };
//!@}

//!@name Static functions
//!@{
  static BondImportance bondImportance(
    const BondIndex& bondIndex,
    const Molecule& molecule,
    const std::map<AtomIndex, unsigned>& smallestCycleMap
  );

  static unsigned distance(const DecisionList& a, const DecisionList& b, const DecisionList& bounds);
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
  //! Generate a new list of discrete dihedral arrangement choices
  DecisionList generateNewDecisionList();
//!@}

//!@name Information
//!@{
  //! Accessor for list of relevant bonds
  const BondList& bondList() const;

  //! Number of conformers needed for full ensemble
  unsigned conformerCount() const;

  //! Try to generate a conformer for
  outcome::result<Utils::PositionCollection> generateConformer(const DecisionList& decisionList);

  //! Infer the decision list from a position collection
  DecisionList getDecisionList(const Utils::PositionCollection& positions);
//!@}

private:
  class Impl;
  std::unique_ptr<Impl> _pImpl;
};

} // namespace molassembler
} // namespace Scine

#endif
