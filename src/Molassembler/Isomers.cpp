/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Isomers.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/Molecule/AtomEnvironmentHash.h"
#include "Molassembler/Graph.h"
#include "Molassembler/AtomStereopermutator.h"
#include "Molassembler/BondStereopermutator.h"
#include "Molassembler/StereopermutatorList.h"
#include "Molassembler/Stereopermutators/AbstractPermutations.h"
#include "Molassembler/Stereopermutators/FeasiblePermutations.h"
#include "Molassembler/Shapes/Data.h"

namespace Scine {
namespace Molassembler {
namespace {

bool partiallyCanonicalizedEnantiomeric(const Molecule& a, const Molecule& b) {
  // Check the precondition
  if(a.graph().N() != b.graph().N() || a.graph().B() != b.graph().B()) {
    return false;
  }

  // Use weaker definition
  constexpr auto bitmask = AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Shapes;

  if(
    a.canonicalComponents() != bitmask
    || b.canonicalComponents() != bitmask
  ) {
    throw std::logic_error(
      "One or more arguments are not partially canonical as required by "
      "the precondition"
    );
  }

  const auto aHashes = Hashes::generate(a.graph().inner(), a.stereopermutators(), bitmask);
  const auto bHashes = Hashes::generate(b.graph().inner(), b.stereopermutators(), bitmask);

  if(aHashes != bHashes) {
    throw std::logic_error("Both molecules are canonical but weaker hashes do not match.");
  }

  /* The number of AtomStereopermutators and BondStereopermutarors must match
   * between both molecules.
   */
  if(
    a.stereopermutators().A() != b.stereopermutators().A()
    && a.stereopermutators().B() != b.stereopermutators().B()
  ) {
    return false;
  }

  /* Comparing stereopermutations for enantiomerism can yield three results:
   * - In this shape, there are no enantiomers
   * - In this shape, there are enantiomers
   *   - These two are enantiomeric
   *   - These two are not enantiomeric
   *
   * If we encounter only pairs in whose shapes no enantiomers exist, then
   * the molecules are not enantiomeric.
   *
   * If we encounter a single pair that is not enantiomeric, then the molecules
   * are not enantiomeric.
   *
   * For the molecules to be enantiomeric, we must encounter at least one
   * enantiomeric pair of permutators and no non-enantiomeric pairs. Pairs
   * of permutators whose shapes have no enantiomers are ignored.
   */
  bool atLeastOneEnantiomericPair = false;

  // Match B's stereopermutators to A's (lists are same-size as established above)
  for(
    const AtomStereopermutator& aPermutator :
    a.stereopermutators().atomStereopermutators()
  ) {
    // Skip all stereopermutators with just one stereopermutation
    if(aPermutator.numStereopermutations() <= 1) {
      continue;
    }

    const auto bPermutatorOption = b.stereopermutators().option(
      aPermutator.placement()
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

    // If the shapes do not match, these molecules cannot be enantiomers
    if(aPermutator.getShape() != bPermutator.getShape()) {
      return false;
    }

    const auto& aPermutation = aPermutator.getAbstract().permutations.list.at(
      *aPermutator.indexOfPermutation()
    );
    const auto& bPermutation = bPermutator.getAbstract().permutations.list.at(
      *bPermutator.indexOfPermutation()
    );

    // Stereopermutation enantiomerism comparison
    boost::optional<bool> enantiomericPairOption = Stereopermutations::enantiomer(
      aPermutation,
      bPermutation,
      aPermutator.getShape()
    );

    if(enantiomericPairOption) {
      if(*enantiomericPairOption) {
        atLeastOneEnantiomericPair = true;
      } else {
        return false;
      }
    }
  }

  if(!atLeastOneEnantiomericPair) {
    return false;
  }

  /* Assignments of bond stereopermutators must match as these are
   * unaffected by mirroring.
   */
  for(
    const BondStereopermutator& permutator :
    a.stereopermutators().bondStereopermutators()
  ) {
    if(permutator.numStereopermutations() <= 1) {
      continue;
    }

    const auto matchOption = b.stereopermutators().option(permutator.placement());

    // There must be a matching permutator
    if(!matchOption) {
      return false;
    }

    // Index of permutation must match
    if(permutator.indexOfPermutation() != matchOption->indexOfPermutation()) {
      return false;
    }
  }

  return true;
}

} // namespace

bool enantiomeric(const Molecule& a, const Molecule& b) {
  constexpr auto componentsBitmask = AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Shapes;

  const auto maybeCanonicalize = [&](const Molecule& m) -> boost::optional<Molecule> {
    if(m.canonicalComponents() != componentsBitmask) {
      boost::optional<Molecule> canonicalCopy = m;
      canonicalCopy->canonicalize(componentsBitmask);
      return canonicalCopy;
    }

    return boost::none;
  };

  const auto maybeCanonicalA = maybeCanonicalize(a);
  const auto maybeCanonicalB = maybeCanonicalize(b);

  return partiallyCanonicalizedEnantiomeric(
    maybeCanonicalA.value_or(a),
    maybeCanonicalB.value_or(b)
  );
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

    // If there are no enantiomers for this shape, we skip the permutator
    const auto& mirrorPermutation = Shapes::mirror(permutator.getShape());
    if(mirrorPermutation.empty()) {
      continue;
    }

    // Find the current permutation
    const auto& currentStereopermutation = permutator.getAbstract().permutations.list.at(
      *permutator.indexOfPermutation()
    );

    // Apply the mirror
    auto mirrored = currentStereopermutation.applyPermutation(mirrorPermutation);

    // Find an existing permutation that is superposable with the mirror permutation
    const auto& permutationsList = permutator.getAbstract().permutations.list;
    auto matchingPermutationIter = std::find_if(
      std::begin(permutationsList),
      std::end(permutationsList),
      [&](const auto& permutation) -> bool {
        return Stereopermutations::rotationallySuperimposable(
          permutation,
          mirrored,
          permutator.getShape()
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

    const auto& feasiblePermutations = permutator.getFeasible().indices;

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
    std::move(stereopermutators)
  };
}

} // namespace Molassembler
} // namespace Scine
