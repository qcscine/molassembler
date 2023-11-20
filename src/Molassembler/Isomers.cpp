/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
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
#include "Molassembler/Temple/Optionals.h"
#include "Molassembler/Temple/Functional.h"

namespace Scine {
namespace Molassembler {
namespace {

bool everythingBesidesStereopermutationsSame(const Molecule& a, const Molecule& b) {
  // Use weaker definition
  constexpr auto bitmask = AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Shapes;

  // Check the precondition
  assert(a.canonicalComponents() == bitmask && b.canonicalComponents() == bitmask);

  // Graph basic equality
  if(a.graph().V() != b.graph().V() || a.graph().E() != b.graph().E()) {
    return false;
  }

  const auto aHashes = Hashes::generate(a.graph().inner(), a.stereopermutators(), bitmask);
  const auto bHashes = Hashes::generate(b.graph().inner(), b.stereopermutators(), bitmask);

  if(aHashes != bHashes) {
    return false;
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

  return true;
}

bool partiallyCanonicalizedEnantiomeric(const Molecule& a, const Molecule& b) {
  if(!everythingBesidesStereopermutationsSame(a, b)) {
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
  return Temple::all_of(
    a.stereopermutators().bondStereopermutators(),
    [&](const BondStereopermutator& permutator) -> bool {
      if(permutator.numStereopermutations() <= 1) {
        return true;
      }

      const auto matchOption = b.stereopermutators().option(permutator.placement());

      if(!matchOption) {
        return false;
      }

      if(permutator.indexOfPermutation() != matchOption->indexOfPermutation()) {
        return false;
      }

      return true;
    }
  );
}

boost::optional<unsigned> permutationDifferences(const Molecule& a, const Molecule &b) {
  unsigned permutationDifferences = 0;
  for(
    const AtomStereopermutator& permutator :
    a.stereopermutators().atomStereopermutators()
  ) {
    auto matchOption = b.stereopermutators().option(permutator.placement());
    if(!matchOption) {
      return boost::none;
    }

    if(permutator.getShape() != matchOption->getShape()) {
      return boost::none;
    }

    if(permutator.indexOfPermutation() != matchOption->indexOfPermutation()) {
      ++permutationDifferences;
    }
  }

  for(
    const BondStereopermutator& permutator :
    a.stereopermutators().bondStereopermutators()
  ) {
    auto matchOption = b.stereopermutators().option(permutator.placement());
    if(!matchOption) {
      return boost::none;
    }

    if(permutator.indexOfPermutation() != matchOption->indexOfPermutation()) {
      ++permutationDifferences;
    }
  }

  return permutationDifferences;
}

bool partiallyCanonicalizedDiastereomeric(const Molecule& a, const Molecule& b) {
  // Different at one or more stereopermutators, but not mirror images
  if(!everythingBesidesStereopermutationsSame(a, b)) {
    return false;
  }

  if(partiallyCanonicalizedEnantiomeric(a, b)) {
    return false;
  }

  return Temple::Optionals::map(
    permutationDifferences(a, b),
    [](unsigned differences) -> bool {
      return differences > 0;
    }
  ).value_or(false);
}

bool partiallyCanonicalizedEpimeric(const Molecule& a, const Molecule& b) {
  // Different at one or more stereopermutators, but not mirror images
  if(!everythingBesidesStereopermutationsSame(a, b)) {
    return false;
  }

  return Temple::Optionals::map(
    permutationDifferences(a, b),
    [](unsigned differences) -> bool {
      return differences == 1;
    }
  ).value_or(false);
}

boost::optional<Molecule> maybeCanonicalize(const Molecule& m) {
  constexpr auto componentsBitmask = AtomEnvironmentComponents::ElementTypes
    | AtomEnvironmentComponents::BondOrders
    | AtomEnvironmentComponents::Shapes;

  if(m.canonicalComponents() != componentsBitmask) {
    boost::optional<Molecule> canonicalCopy = m;
    canonicalCopy->canonicalize(componentsBitmask);
    return canonicalCopy;
  }

  return boost::none;
}

} // namespace

bool enantiomeric(const Molecule& a, const Molecule& b) {
  const auto maybeCanonicalA = maybeCanonicalize(a);
  const auto maybeCanonicalB = maybeCanonicalize(b);

  return partiallyCanonicalizedEnantiomeric(
    maybeCanonicalA.value_or(a),
    maybeCanonicalB.value_or(b)
  );
}

bool diastereomeric(const Molecule& a, const Molecule& b) {
  const auto maybeCanonicalA = maybeCanonicalize(a);
  const auto maybeCanonicalB = maybeCanonicalize(b);

  return partiallyCanonicalizedDiastereomeric(
    maybeCanonicalA.value_or(a),
    maybeCanonicalB.value_or(b)
  );
}

bool epimeric(const Molecule& a, const Molecule& b) {
  const auto maybeCanonicalA = maybeCanonicalize(a);
  const auto maybeCanonicalB = maybeCanonicalize(b);

  return partiallyCanonicalizedEpimeric(
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
    const auto mirrored = currentStereopermutation.applyPermutation(mirrorPermutation);

    const auto matchingAssignment = Stereopermutators::Feasible::findRotationallySuperposableAssignment(
      currentStereopermutation,
      permutator.getShape(),
      permutator.getAbstract(),
      permutator.getFeasible()
    );

    // NOTE: This unassigns the stereopermutator if no matching assignment was found
    permutator.assign(matchingAssignment);
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
