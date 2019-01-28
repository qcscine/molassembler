/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "Isomers.h"

#include "molassembler/Molecule.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Stereopermutators/PermutationState.h"

namespace Scine {

namespace molassembler {

bool enantiomeric(const Molecule& /* a */, const Molecule& /* b */) {
  /* Two molecules are mirror images of one another if a graph isomorphism
   * that matches
   * - element type
   * - bond type
   * - symmetry
   * - stereopermutations (of BondStereopermutators only!)
   *
   * yields a mapping where all pairs of atom stereopermutators have
   * mirror assignments.
   */

  throw std::runtime_error("Not implemented!");
}

bool enantiomeric(
  const Molecule& a,
  const Molecule& b,
  SameIndexingTag /* sameIndexingTag */
) {
  // Check the precondition
  if(a.graph().N() != b.graph().N()) {
    throw std::logic_error(
      "Enantiomer check precondition violated: Molecules are not of same size."
    );
  }
  constexpr auto bitmask = AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Symmetries;

  const auto aHashes = hashes::generate(a.graph().inner(), a.stereopermutators(), bitmask);
  const auto bHashes = hashes::generate(b.graph().inner(), b.stereopermutators(), bitmask);

  if(aHashes != bHashes) {
    throw std::logic_error(
      "Enantiomer check precondition violated: Molecules are not "
      "indexed identically."
    );
  }

  /* The number of AtomStereopermutators must match between both molecules.
   *
   * We don't care about BondStereopermutators since these are unchanged by
   * mirroring operations.
   *
   * TODO that isn't strictly true as soon as bond stereopermutators are
   * more complicated than 2x trigonal planar
   */
  if(a.stereopermutators().A() != b.stereopermutators().A()) {
    return false;
  }

  /* Comparing stereopermutations for enantiomerism can yield three results:
   * - In this symmetry, there are no enantiomers
   * - In this symmetry, there are enantiomers
   *   - These two are enantiomeric
   *   - These two are not enantiomeric
   *
   * If we encounter only pairs in whose symmetries no enantiomers exist, then
   * the molecules are not enantiomeric.
   *
   * If we encounter a single pair that is not enantiomeric, then the molecules
   * are not enantiomeric.
   *
   * For the molecules to be enantiomeric, we must encounter at least one
   * enantiomeric pair of permutators and no non-enantiomeric pairs. Pairs
   * of permutators whose symmetries have no enantiomers are ignored.
   */
  bool atLeastOneEnantiomericPair = false;

  // Match B's stereopermutators to A's (lists are same-size as established above)
  for(const AtomStereopermutator& aPermutator : a.stereopermutators().atomStereopermutators()) {
    // Skip all stereopermutators with just one stereopermutation
    if(aPermutator.numStereopermutations() <= 1) {
      continue;
    }

    const auto& bPermutatorOption = b.stereopermutators().option(
      aPermutator.centralIndex()
    );

    // If there is no matching permutator, then these cannot be enantiomers
    if(!bPermutatorOption) {
      return false;
    }

    const AtomStereopermutator& bPermutator = *bPermutatorOption;

    /* If one or both permutators are unassigned, then we cannot say for certain
     * if a 3D representation of both will be enantiomeric, hence we return
     * false
     */
    if(!aPermutator.assigned() || !bPermutator.assigned()) {
      return false;
    }

    // If the symmetries do not match, these molecules cannot be enantiomers
    if(aPermutator.getSymmetry() != bPermutator.getSymmetry()) {
      return false;
    }

    const auto& aPermutation = aPermutator.getPermutationState().permutations.assignments.at(
      *aPermutator.indexOfPermutation()
    );
    const auto& bPermutation = bPermutator.getPermutationState().permutations.assignments.at(
      *bPermutator.indexOfPermutation()
    );

    // Stereopermutation enantiomerism comparison
    boost::optional<bool> enantiomericPairOption = aPermutation.isEnantiomer(
      bPermutation,
      aPermutator.getSymmetry()
    );

    if(enantiomericPairOption) {
      if(*enantiomericPairOption) {
        atLeastOneEnantiomericPair = true;
      } else {
        return false;
      }
    }
  }

  return atLeastOneEnantiomericPair;
}

Molecule enantiomer(const Molecule& a) {
  /* Copy a's StereopermutatorList, then find and set the enantiomeric
   * stereopermutation for each atom stereopermutator with more than one
   * stereopermutation
   */
  StereopermutatorList stereopermutators = a.stereopermutators();

  for(auto& permutator : stereopermutators.atomStereopermutators()) {
    /* If there is only a single permutation or the permutator is unassigned,
     * we leave it be
     */
    if(permutator.numStereopermutations() <= 1 || !permutator.assigned()) {
      continue;
    }

    // If there are no enantiomers for this symmetry, we skip the permutator
    const auto& mirrorPermutation = Symmetry::mirror(permutator.getSymmetry());
    if(mirrorPermutation.empty()) {
      continue;
    }

    // Find the current permutation
    auto currentStereopermutation = permutator.getPermutationState().permutations.assignments.at(
      *permutator.indexOfPermutation()
    );

    // Apply the mirror
    currentStereopermutation.applyRotation(mirrorPermutation);

    // Find an existing permutation that is superposable with the mirror permutation
    const auto& permutationsList = permutator.getPermutationState().permutations.assignments;
    auto matchingPermutationIter = std::find_if(
      std::begin(permutationsList),
      std::end(permutationsList),
      [&](const auto& permutation) -> bool {
        return permutation.isRotationallySuperimposable(
          currentStereopermutation,
          permutator.getSymmetry()
        );
      }
    );

    // If we cannot find a matching stereopermutation, then unassign
    if(matchingPermutationIter == std::end(permutationsList)) {
      permutator.assign(boost::none);
      continue;
    }

    // Now we have a stereopermutation index, but we need an assignment index
    unsigned stereopermutationIndex = matchingPermutationIter - std::begin(permutationsList);

    const auto& feasiblePermutations = permutator.getPermutationState().feasiblePermutations;

    auto assignmentIter = std::find(
      std::begin(feasiblePermutations),
      std::end(feasiblePermutations),
      stereopermutationIndex
    );

    // If the found permutation is infeasible, then unassign
    if(assignmentIter == std::end(feasiblePermutations)) {
      permutator.assign(boost::none);
      continue;
    }

    // Otherwise, we can assign it with the found permutation
    permutator.assign(assignmentIter - std::begin(feasiblePermutations));
  }

  /* Construct a new molecule by copying the graph and moving in the altered
   * permutator list
   */
  return {
    a.graph(),
    std::move(stereopermutators),
    AtomEnvironmentComponents::None
  };
}

} // namespace molassembler

} // namespace Scine
