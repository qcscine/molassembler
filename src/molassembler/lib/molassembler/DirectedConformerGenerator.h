/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Directed conformer generation class and helper functions
 */

#ifndef INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H
#define INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H

#include "molassembler/Conformers.h"

#include "boost/optional/optional_fwd.hpp"
#include "boost/variant/variant_fwd.hpp"
#include <map>
#include <memory>

namespace Scine {
namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

class Molecule;
class BondStereopermutator;

/**
 * @brief Helper type for directed conformer generation.
 *
 * Generates guaranteed new combinations of BondStereopermutator assignments
 * and provides helper functions for the generation of conformers using these
 * combinations and the reverse, finding the combinations from conformers.
 *
 * @note This type is not copyable.
 */
class DirectedConformerGenerator {
public:
//!@name Public types
//!@{
  //! Type used to represent the list of bonds relevant to directed conformer generation
  using BondList = std::vector<BondIndex>;
  /*!
   * @brief Type used to represent assignments at bonds
   * @note You can serialize / deserialize this with Scine::base64::encode and
   *   Scine::base64::decode. It's not the most efficient representation
   *   but still better than each position having its own character.
   */
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
   * The distance metric is:
   *
   * \f$d = \sum_i \min\left(a_i - b_i \textrm{ mod } U_i, b_i - a_i \textrm{ mod } U_i\right)\f$
   *
   * where \f$a_i\f$ is the choice in @p a at position \f$i\f$ and likewise for
   * @p b, and \f$U_i\f$ is the upper exclusive bound on the choice values at
   * position \f$i\f$.
   *
   * This is akin to the shortest distance between the choices when arranged in
   * a modular number circle.
   *
   * \code{cpp}
   * auto a = std::vector<std::uint8_t> {{1, 5}};
   * auto b = std::vector<std::uint8_t> {{3, 0}};
   * auto bounds = std::vector<std::uint8_t> {{5, 6}};
   * unsigned d = distance(a, b, bounds); // Yields 2 + 1 = 3
   * \endcode
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
  /**
   * @brief Constructor
   *
   * @param molecule Molecule for which to generate conformers
   * @param bondsToConsider A list of suggestions of which bonds to consider.
   *   Bonds for which considerBond yields an IgnoreReason will still be
   *   ignored. If the list is empty, all bonds of a molecule will be
   *   considered.
   *
   * Scales linearly with the number of bonds in @p molecule or
   * @p bondsToConsider's size.
   */
  explicit DirectedConformerGenerator(
    Molecule molecule,
    const BondList& bondsToConsider = {}
  );
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
   *
   * Guarantees that the generated list is not yet part of the underlying set.
   * Scales linearly with the number of considered dihedrals.
   *
   * @throws std::logic_error If the underlying set is full, i.e. all decision
   *   lists for conformers have been generated.
   * @post The new DecisionList is part of the stored list of generated
   *   decision lists and will not be generated again.
   */
  DecisionList generateNewDecisionList();

  /*!
   * @brief Adds a decision list to the underlying set-like data structure
   *
   * Scales linearly with the length of @p decisionList.
   *
   * @returns @p true if @p decisionList wasn't already part of the set
   */
  bool insert(const DecisionList& decisionList);

  /*!
   * @brief Checks whether a DecisionList is part of the underlying set
   *
   * Scales linearly with the length of @p decisionList.
   */
  bool contains(const DecisionList& decisionList);
//!@}

//!@name Information
//!@{
  //! Accessor for list of relevant bonds, O(1)
  const BondList& bondList() const;

  //! Number of conformer decision lists stored in the underlying set-like data structure, O(1)
  unsigned conformerCount() const;

  //! Number of conformers needed for full ensemble, O(1)
  unsigned idealEnsembleSize() const;

  /*!
   * @brief Try to generate a conformer for a particular decision list
   *
   * This is very similar to the free generateConformation function in terms
   * of what @p configuration will accept.
   */
  outcome::result<Utils::PositionCollection> generateConformation(
    const DecisionList& decisionList,
    const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
  );

  /*!
   * @brief Infer a decision list from positional information
   *
   * @warning This function assumes several things about your supplied positions
   * - There have only been dihedral changes and no AtomStereopermutator
   *   assignment changes
   * - The molecule represented in @p positions has not constutitionally
   *   rearranged
   *
   * @throws std::logic_error If an assignment could not be recovered from
   *   positions
   */
  DecisionList getDecisionList(Utils::PositionCollection positions) const;
//!@}

private:
  class Impl;
  std::unique_ptr<Impl> _pImpl;
};

} // namespace molassembler
} // namespace Scine

#endif
