/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Directed conformer generation class and helper functions
 */

#ifndef INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H
#define INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H

#include "molassembler/Conformers.h"
#include "molassembler/BondStereopermutator.h"

#include "boost/optional/optional_fwd.hpp"
#include "boost/variant/variant_fwd.hpp"
#include "boost/functional/hash.hpp"

#include <unordered_map>
#include <memory>

namespace Scine {

namespace Utils {
class AtomCollection;
} // namespace Utils

namespace molassembler {

namespace outcome = BOOST_OUTCOME_V2_NAMESPACE;

class Molecule;
class BondStereopermutator;

/** @brief Helper type for directed conformer generation.
 *
 * Generates guaranteed new combinations of BondStereopermutator assignments
 * and provides helper functions for the generation of conformers using these
 * combinations and the reverse, finding the combinations from conformers.
 *
 * Client code could look something like this. Make sure to differentiate
 * cases in which the list of considered bonds by bondList() is empty, since
 * many member functions behave differently in those circumstances.
 *
 * @code{cpp}
 * auto mol = IO::read(...);
 * std::vector<Utils::PositionCollection> conformers;
 * DirectedConformerGenerator generator {mol};
 * if(generator.bondList().empty()) {
 *   // The generator has decided that there are no bonds that need to be
 *   // systematically rotated for directed conformer generation. You might
 *   // as well just use a conformation directly from generateConformation for
 *   // the sole conformation needed for a full ensemble
 *   auto conformerResult = generateConformation(mol);
 *   if(conformerResult) {
 *     conformers.push_back(std::move(conformerResult.value()));
 *   } else {
 *     std::cout << "Could not generate conformer: " << conformerResult.error().message() << "\n";
 *   }
 * } else {
 *   while(generator.decisionListSetSize() != generator.idealEnsembleSetSize()) {
 *     auto newDecisionList = generator.generateNewDecisionList();
 *     auto conformerResult = generator.generateConformation(newDecisionList);
 *     if(conformerResult) {
 *       conformers.push_back(std::move(conformerResult.value()));
 *     } else {
 *       std::cout << "Could not generate conformer: " << conformerResult.error().message() << "\n";
 *     }
 *   }
 * }
 * @endcode
 *
 * @note This type is not copyable.
 */
class DirectedConformerGenerator {
public:
//!@name Public types
//!@{
  //! Type used to represent the list of bonds relevant to directed conformer generation
  using BondList = std::vector<BondIndex>;
  /*! @brief Type used to represent assignments at bonds
   *
   * @note You can serialize / deserialize this with Scine::base64::encode and
   *   Scine::base64::decode. It's not the most efficient representation
   *   but still better than each position having its own character.
   */
  using DecisionList = std::vector<std::uint8_t>;

  //! Value set in decision lists if no decision could be recovered
  constexpr static std::uint8_t unknownDecision = std::numeric_limits<std::uint8_t>::max();

  //* @brief Reason why a bond is ignored
  enum class IgnoreReason {
    //! There is not an assigned stereopermutator on both ends of the bond
    AtomStereopermutatorPreconditionsUnmet,
    //! There is already an assigned bond stereopermutator on the bond
    HasAssignedBondStereopermutator,
    //! At least one constituting atom is terminal
    HasTerminalConstitutingAtom,
    /*! @brief This bond is in a cycle
     *
     * Despite the fact that cycle bonds may very well contribute to the
     * conformational ensemble, it is difficult to reason about conformational
     * flexibility of cycles:
     * - Is a cycle aromatic or anti-aromatic?
     * - Is there a partial conjugated system?
     * - Are trans-arrangements of cycle atom sequences feasible in the cycle?
     *
     * We see four possible strategies of dealing with cycle bonds:
     * 1. Consider all of them.
     * 2. Add chemical intuitive reasoning to exclude some bonds in cycles from
     *    consideration
     * 3. Use Distance Geometry to reason about possible dihedrals
     * 4. Consider none of them.
     *
     * The first strategy has several important drawbacks: Dihedrals are heavily
     * restricted and/or correlated in the chemically common small cycles, and
     * most, if not nearly all, combinations will not be representable in three
     * dimensions. Bonds from cycles incur heavy cost in the decision list set
     * representation and computational time needed to generate conformations
     * because it is not possible with mere triangle inequality smoothing to
     * determine representability of these conformers, and hence a full
     * refinement is done for each.
     *
     * In contrast, the second strategy could yield properly limited assignment
     * possibilities for common chemical patterns. The algorithms needed to
     * answer the questions listed above are complex and could easily fail
     * outside common organic chemical patterns.
     *
     * The last strategy is cleanest, but also likely considerably
     * computationally expensive.
     *
     * For now, we take strategy number four - ignoring bonds in cycles for
     * directed conformer generation - until we can dedicate some resources to
     * approach three.
     */
    InCycle,
    /*! @brief This bond is an eta bond (indicates bonding to haptic ligands,
     *  and therefore excluded)
     */
    IsEtaBond,
    /*! @brief Rotation about this bond is isotropic (all ligands have same
     *  ranking on at least one side)
     */
    RotationIsIsotropic
  };
//!@}

//!@name Static functions
//!@{
  /** @brief Decide whether to consider a bond's dihedral values for directed
   *   conformer generation or not
   *
   * @param bondIndex The bond to consider
   * @param molecule The molecule in which @p bond exists
   * @param smallestCycleMap A map of atom indices to the smallest cycle
   *   they are in
   *
   * @complexity{@math{O(S!)} where @math{S} is the size of the larger shape
   * constituting @p bondIndex}
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
    const std::unordered_map<AtomIndex, unsigned>& smallestCycleMap
  );

  /** @brief Calculates a distance metric between two decision lists for
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
   * @complexity{@math{\Theta(N)}}
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
  /** @brief Constructor
   *
   * @param molecule Molecule for which to generate conformers
   * @param bondsToConsider A list of suggestions of which bonds to consider.
   *   Bonds for which considerBond() yields an IgnoreReason will still be
   *   ignored. If the list is empty, all bonds of a molecule will be
   *   tested against considerBond().
   *
   * @complexity{@math{\Theta(B)} where @math{B} is the number of bonds in the
   * molecule. If there is a particularly large shape in the molecule, this
   * can dominate with @math{\Theta(S!)}.}
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
  /*! @brief Generate a new list of discrete dihedral arrangement choices
   *
   * Guarantees that the generated list is not yet part of the underlying set.
   *
   * @complexity{@math{\Theta(N)}}
   *
   * @throws std::logic_error If the underlying set is full, i.e. all decision
   *   lists for conformers have been generated or if there are no bonds to
   *   consider.
   * @post The new DecisionList is part of the stored list of generated
   *   decision lists and will not be generated again. The result of
   *   decisionListSetSize() is incremented.
   *
   * @returns a DecisionList of length matching the number of relevant bonds.
   */
  DecisionList generateNewDecisionList();

  /*!
   * @brief Adds a decision list to the underlying set-like data structure
   *
   * @complexity{@math{\Theta(N)}}
   *
   * @throws std::logic_error If the result of bondList() is empty, i.e. there
   *   are no bonds to consider for directed conformer generation.
   *
   * @returns @p true if @p decisionList wasn't already part of the set
   */
  bool insert(const DecisionList& decisionList);

  /*!
   * @brief Checks whether a DecisionList is part of the underlying set
   *
   * @throws std::logic_error If the result of bondList() is empty, i.e. there
   *   are no bonds to consider for directed conformer generation.
   *
   * @complexity{@math{\Theta(N)}}
   */
  bool contains(const DecisionList& decisionList);
//!@}

//!@name Information
//!@{
  /*!
   * @brief Accessor for list of relevant bonds
   *
   * @complexity{@math{\Theta(1)}}
   * @note This list may be empty. Many member functions may throw under these
   *   conditions.
   */
  const BondList& bondList() const;

  /*!
   * @brief Number of conformer decision lists stored in the underlying
   *   set-like data structure
   * @returns The number of DecisionLists stored in the underlying set.
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @warning If bondList() returns an empty list, i.e. there are no bonds to
   *   consider for directed conformer generation, this always returns zero.
   */
  unsigned decisionListSetSize() const;

  /*! @brief Number of conformers needed for full ensemble
   *
   * @complexity{@math{\Theta(1)}}
   *
   * @warning If bondList() returns an empty list, i.e. there are no bonds to
   *   consider for directed conformer generation, this always returns zero.
   */
  unsigned idealEnsembleSize() const;

  /*! @brief Try to generate a conformer for a particular decision list
   *
   * This is very similar to the free generateRandomConformation function in
   * terms of what @p configuration will accept.
   *
   * @see Scine::molassembler::generateRandomConformation()
   *
   * @note Advances the state of the global PRNG.
   *
   * @throws std::invalid_argument If the passed decisionList does not match
   *   the length of the result of bondList().
   */
  outcome::result<Utils::PositionCollection> generateRandomConformation(
    const DecisionList& decisionList,
    const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
  );

  /*! @brief Try to generate a conformer for a particular decision list
   *
   * This is very similar to the free generateConformation function in terms
   * of what @p configuration will accept.
   *
   * @see Scine::molassembler::generateConformation()
   *
   * @throws std::invalid_argument If the passed decisionList does not match
   *   the length of the result of bondList().
   */
  outcome::result<Utils::PositionCollection> generateConformation(
    const DecisionList& decisionList,
    const unsigned seed,
    const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {}
  );

  /*! @brief Yields a molecule reference for a particular decision list
   *
   * @complexity{@math{\Theta(N)} bond stereopermutator assignments}
   */
  const Molecule& conformationMolecule(const DecisionList& decisionList);

  /*! @brief Infer a decision list for relevant bonds from an atom collection
   *
   * For all bonds considered relevant (i.e. all bonds in bondList()), fits
   * supplied positions to possible stereopermutations and returns the result.
   * Entries have the value `DirectedConformer::unknownDecision` if no
   * permutation could be recovered. The usual BondStereopermutator fitting
   * tolerances apply.
   *
   * @warning This function assumes several things about your supplied positions
   * - There have only been dihedral changes and no AtomStereopermutator
   *   assignment changes
   * - The molecule represented in @p positions has not constutitionally
   *   rearranged (although a little check for matching element types does
   *   exist here. This is not a full safeguard against index permutations.)
   *
   * @complexity{@math{\Theta(N)} bond stereopermutator fits}
   *
   * @throws std::logic_error If the element type sequence of the atom
   * collection does not match the underlying molecule
   */
  DecisionList getDecisionList(
    const Utils::AtomCollection& atomCollection,
    BondStereopermutator::FittingMode mode = BondStereopermutator::FittingMode::Thresholded
  );

  /*! @brief Infer a decision list for relevant bonds from positional information only
   *
   * For all bonds considered relevant (i.e. all bonds in bondList()), fits
   * supplied positions to possible stereopermutations and returns the result.
   * Entries have the value `DirectedConformer::unknownDecision` if no
   * permutation could be recovered. The usual BondStereopermutator fitting
   * tolerances apply.
   *
   * @warning This function assumes several things about your supplied positions
   * - There have only been dihedral changes and no AtomStereopermutator
   *   assignment changes
   * - The molecule represented in @p positions has not constutitionally
   *   rearranged
   *
   * @complexity{@math{\Theta(N)} bond stereopermutator fits}
   */
  DecisionList getDecisionList(
    const Utils::PositionCollection& positions,
    BondStereopermutator::FittingMode mode = BondStereopermutator::FittingMode::Thresholded
  );
//!@}

private:
  class Impl;
  std::unique_ptr<Impl> _pImpl;
};

} // namespace molassembler
} // namespace Scine

#endif
