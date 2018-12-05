#include "Isomers.h"

#include "molassembler/Molecule.h"
#include "molassembler/Molecule/AtomEnvironmentHash.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/StereopermutatorList.h"
#include "molassembler/Stereopermutators/PermutationState.h"

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
  constexpr auto bitmask = temple::make_bitmask(AtomEnvironmentComponents::ElementTypes)
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
   * TODO this isn't strictly true as soon as bond stereopermutators are
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

} // namespace molassembler
