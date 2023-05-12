/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Directed conformer generation class and helper functions
 */

#ifndef INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H
#define INCLUDE_MOLASSEMBLER_DIRECTED_CONFORMER_GENERATOR_H

#include "Molassembler/Conformers.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/Options.h"

#include "boost/variant/variant_fwd.hpp"

#include <unordered_map>

namespace Scine {
namespace Utils {
class AtomCollection;
} // namespace Utils

namespace Molassembler {

class Molecule;

/** @brief Helper type for directed conformer generation.
 *
 * Generates new combinations of BondStereopermutator assignments and provides
 * helper functions for the generation of conformers using these combinations
 * and the reverse, finding the combinations from conformers.
 *
 * It is important that you lower your expectations for the modeling of
 * dihedral energy minima, however. Considering that Molassembler neither
 * requires you to supply a correct graph, never detects or kekulizes aromatic
 * systems nor asks you to supply an overall charge for a molecule, it should
 * be understandable that the manner in which Molassembler decides where
 * dihedral energy minima are is somewhat underpowered. The manner in which
 * shape vertices are aligned in stereopermutation enumeration isn't even
 * strictly based on a physical principle. We suggest the following to make the
 * most of what the library can do for you:
 *
 * - Read the documentation for the various alignments. Consider using not just
 *   the default Staggered alignment, but either EclipsedAndStaggered or
 *   BetweenEclipsedAndStaggered to improve your chances of capturing all
 *   rotational minima. This will likely generate more conformers than strictly
 *   required, but should capture all minima.
 * - Energy minimize all generated conformers with a suitable method and then
 *   deduplicate.
 * - Consider using the Relabeler to do a final deduplication step.
 *
 * Client code for generating conformers could look something like this. Make
 * sure to differentiate cases in which the list of considered bonds by
 * bondList() is empty, since many member functions behave differently in those
 * circumstances.
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
class MASM_EXPORT DirectedConformerGenerator {
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

  struct Relabeler;
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
   * @param alignment Alignment to generate BondStereopermutator instances with.
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
    BondStereopermutator::Alignment alignment = BondStereopermutator::Alignment::Staggered
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
   * @param alignment Alignment with which to generate BondStereopermutator on
   *   considered bonds
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
    BondStereopermutator::Alignment alignment = BondStereopermutator::Alignment::Staggered,
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
   * @parblock @note This function advances the state of the global PRNG if the
   *   default argument for @p engine is chosen.
   * @endparblock
   *
   * @returns a DecisionList of length matching the number of relevant bonds.
   */
  DecisionList generateNewDecisionList(Random::Engine& engine = randomnessEngine());

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
  bool contains(const DecisionList& decisionList) const;
//!@}

//!@name Information
//!@{
  //! Get alignment with which this generator was instantiated with
  BondStereopermutator::Alignment alignment() const;

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
   * @see Scine::Molassembler::generateRandomConformation()
   *
   * @note Advances the state of the global PRNG. Not reentrant.
   *
   * @throws std::invalid_argument If the passed decisionList does not match
   *   the length of the result of bondList().
   */
  Result<Utils::PositionCollection> generateRandomConformation(
    const DecisionList& decisionList,
    const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {},
    BondStereopermutator::FittingMode fitting = BondStereopermutator::FittingMode::Nearest
  ) const;

  /*! @brief Try to generate a conformer for a particular decision list
   *
   * This is very similar to the free generateConformation function in terms
   * of what @p configuration will accept.
   *
   * @see Scine::Molassembler::generateConformation()
   *
   * @throws std::invalid_argument If the passed decisionList does not match
   *   the length of the result of bondList().
   */
  Result<Utils::PositionCollection> generateConformation(
    const DecisionList& decisionList,
    unsigned seed,
    const DistanceGeometry::Configuration& configuration = DistanceGeometry::Configuration {},
    BondStereopermutator::FittingMode fitting = BondStereopermutator::FittingMode::Nearest
  ) const;

  /*! @brief Yields a molecule reference for a particular decision list
   *
   * @complexity{@math{\Theta(N)} bond stereopermutator assignments}
   */
  Molecule conformationMolecule(const DecisionList& decisionList) const;

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
    BondStereopermutator::FittingMode mode = BondStereopermutator::FittingMode::Nearest
  ) const;

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
  ) const;

  //! @brief Settings for enumeration
  struct EnumerationSettings {
    EnumerationSettings() : configuration() {}

    //! Number of attempts to generate the correct decision list
    unsigned dihedralRetries = 3;
    //! How a decision list is fitted to a generated conformer
    BondStereopermutator::FittingMode fitting = BondStereopermutator::FittingMode::Nearest;
    //! Conformer generation settings
    DistanceGeometry::Configuration configuration;
  };

  /*! @brief Enumerate all conformers of the captured molecule
   *
   * Clears the stored set of decision lists, then enumerates all conformers of
   * the molecule in parallel.
   *
   * @param callback Function called with decision list and conformer
   *    positions for each successfully generated. It is guaranteed that the
   *    callback function is never called simultaneously even in parallel
   *    execution.
   * @param seed Randomness initiator for decision list and conformer
   *    generation
   * @param settings Further parameters for enumeration algorithms
   *
   * @parblock @note This function is parallelized. Use the OMP_NUM_THREADS
   * environment variable to control the number of threads used. Callback
   * invocations are unsequenced but the arguments are reproducible.
   * @endparblock
   */
  void enumerate(
    std::function<void(const DecisionList&, Utils::PositionCollection)> callback,
    unsigned seed,
    const EnumerationSettings& settings = {}
  );

  /*! @brief Enumerate all conformers of the captured molecule
   *
   * Clears the stored set of decision lists, then enumerates all conformers of
   * the molecule in parallel.
   *
   * @param callback Function called with decision list and conformer
   *   positions for each successfully generated pair. It is guaranteed that
   *   the callback function is never called simultaneously even in parallel
   *   execution.
   * @param seed Randomness initiator for decision list and conformer
   *   generation
   * @param settings Further parameters for enumeration algorithms
   *
   * @parblock @note This function is parallelized. Use the OMP_NUM_THREADS
   * environment variable to control the number of threads used. Callback
   * invocations are unsequenced but the arguments are reproducible.
   * @endparblock
   *
   * @parblock @note This function advances the state of the global PRNG.
   * @endparblock
   */
  void enumerateRandom(
    std::function<void(const DecisionList&, Utils::PositionCollection)> callback,
    const EnumerationSettings& settings = {}
  );

  //! Generates a relabeler for the molecule and considered bonds
  Relabeler relabeler() const;

  //! Relabels a DecisionList into bin midpoint integers
  //! Returns the dihedra angles. The angle is defined as the first angle of the
  //! stereopermutation at the given bond (note that this is inconsistent
  //! with its definition in Relabeler::add).
  [[deprecated]]
  std::vector<int> binMidpointIntegers(const DecisionList& decision) const;

  //! Relabels a DecisionList into the bounds of its bin
  std::vector<std::pair<int, int>> binBounds(const DecisionList& decision) const;
//!@}

private:
  class Impl;
  std::unique_ptr<Impl> pImpl_;
};

/*! @brief Relabeler for decision lists with minimized structures
 *
 * Type to help with relabeling decision lists with true minima bins
 * once a set of conformers has been energy minimized with a suitable
 * procedure. This helps to deal with dihedral extremum placements that differ
 * significantly from Molassembler's simple expectations.
 *
 * Rotational symmetries from the hypothesized alignments are reused to map
 * rotationally equivalent dihedrals down to single bins.
 */
struct DirectedConformerGenerator::Relabeler {
//!@name Types
//!@{
  struct DihedralInfo {
    std::vector<AtomIndex> is;
    AtomIndex j;
    AtomIndex k;
    std::vector<AtomIndex> ls;
    unsigned symmetryOrder;
  };
  using Interval = std::pair<double, double>;
  using Intervals = std::vector<Interval>;
//!@}

  /*! @brief Simplest density-based binning function
   *
   * Just sorts the dihedral values and then considers any values within the
   * delta as part of the same bin.
   */
  static Intervals densityBins(
    const std::vector<double>& dihedrals,
    double delta,
    unsigned symmetryOrder = 1
  );

  static std::pair<double, double> makeBounds(double phi, double tolerance);
  static std::pair<int, int> integerBounds(const std::pair<double, double>& bounds);

  /*! @brief Construct a relabeler with a custom list of bonds
   *
   * @param bonds List of bonds to consider.
   * @param mol Molecule whose bonds we want to consider. Needs to have
   * BondStereopermutators instantiated on bonds in @p bonds.
   *
   * @note A relabeler with the inferred list of bonds that can be considered
   * as determined by a DirectedConformerGenerator can be obtained by calling
   * DirectedConformerGenerator::relabeler().
   */
  Relabeler(const DirectedConformerGenerator::BondList& bonds, const Molecule& mol);

  //! Add a structure to the set to relabel. Yields the considered dihedrals
  std::vector<double> add(const Utils::PositionCollection& positions);

  /*! Generate bins for each set of observed dihedrals
   *
   * Yields a vector of dimension equal to the bonds, with each element a
   * vector with the intervals for the bond.
   */
  std::vector<Intervals> bins(double delta=M_PI / 6) const;

  /*! @brief Determine relabeling for all added position sets in order
   *
   * Call this as soon as all positions to be reclassified have been added.
   *
   * Returns relabeling for each set of positions in order.
   *
   * Yields a vector of dimension equal to the number of structures, with each
   * element a vector of dimension equal to the number of bonds.
   */
  std::vector<std::vector<unsigned>> binIndices(
    const std::vector<Intervals>& allBins
  ) const;

  //! Relabel bin indices for all structures with bin midpoint integers
  std::vector<std::vector<int>> binMidpointIntegers(
    const std::vector<std::vector<unsigned>>& binIndices,
    const std::vector<Intervals>& allBins
  ) const;

  //! Relabel bin indices for all structures with bin bounds
  std::vector<
    std::vector<std::pair<int, int>>
  > binBounds(
    const std::vector<std::vector<unsigned>>& binIndices,
    const std::vector<Intervals>& allBins
  ) const;

//!@name State
//!@{
  std::vector<DihedralInfo> sequences;
  std::vector<std::vector<double>> observedDihedrals;
//!@}
};

} // namespace Molassembler
} // namespace Scine

#endif
