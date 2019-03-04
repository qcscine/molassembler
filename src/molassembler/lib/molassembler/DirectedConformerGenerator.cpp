#include "DirectedConformerGenerator.h"

#include "molassembler/Cycles.h"
#include "molassembler/Molecule.h"
#include "molassembler/OuterGraph.h"
#include "molassembler/StereopermutatorList.h"

#include "stereopermutation/Composites.h"
#include "temple/BoundedNodeTrie.h"
#include "temple/Functional.h"

namespace Scine {
namespace molassembler {

namespace detail {

int euclideanModulo(const int a, const int base) {
  return ((a % base) + base) % base;
}

unsigned distance(const int i, const int j, const int U) {
  assert(i >= 0 && j >= 0 && U > 1);

  int value = std::min(
    euclideanModulo(i - j, U),
    euclideanModulo(j - i, U)
  );

  assert(value > 0);

  return static_cast<unsigned>(value);
}

} // namespace detail

class DirectedConformerGenerator::Impl {
public:
  Impl(Molecule molecule) : _molecule(std::move(molecule)) {
    // Figure out which bonds are relevant, in order (!)
    // Collect bounds for each
    // initialize _decisionLists with the bounds
  }

private:
  Molecule _molecule;
  BondList _relevantBonds;
  DecisionList _bounds;
  temple::BoundedNodeTrie<std::uint8_t> _decisionLists;
};

/* Static functions */
unsigned DirectedConformerGenerator::distance(const DecisionList& a, const DecisionList& b, const DecisionList& bounds) {
  if(a.size() != b.size() || a.size() != bounds.size()) {
    throw std::invalid_argument("Not all decision lists have the same length");
  }

  unsigned distance = 0;
  for(unsigned i = 0; i < bounds.size(); ++i) {
    assert(a[i] < bounds[i] && b[i] < bounds[i]);
    distance += detail::distance(a[i], b[i], bounds[i]);
  }

  return distance;
}

DirectedConformerGenerator::BondImportance DirectedConformerGenerator::bondImportance(
  const BondIndex& bondIndex,
  const Molecule& molecule,
  const std::map<AtomIndex, unsigned>& smallestCycleMap
) {
  // Make sure the bond exists in the first place
  assert(molecule.graph().adjacent(bondIndex.first, bondIndex.second));

  /* The dihedral already has an assigned bond stereopermutator placed on it or
   * that bond stereopermutator is isotropic (ranking of all substituents at
   * one side is the same)
   */
  if(auto bondStereopermutatorOption = molecule.stereopermutators().option(bondIndex)) {
    if(bondStereopermutatorOption->assigned()) {
      return BondImportance::IgnoreHasAssignedBondStereopermutator;
    }

    if(bondStereopermutatorOption->composite().isIsotropic()) {
      return BondImportance::IgnoreRotationIsIsotropic;
    }
  }

  /* Check with cycles information next: If the dihedral is in a cycle of size 4
   * or lower, we will not consider the dihedral at all
   */
  auto inCycleOfSizeThreeOrFour = [&smallestCycleMap](const AtomIndex i) -> bool {
    auto findIter = smallestCycleMap.find(i);
    if(findIter != std::end(smallestCycleMap)) {
      return findIter->second <= 4;
    }

    return false;
  };

  if(inCycleOfSizeThreeOrFour(bondIndex.first) || inCycleOfSizeThreeOrFour(bondIndex.second)) {
    return BondImportance::IgnoreInSmallCycle;
  }

  /* Last, instantiate a BondStereopermutator on the bond and use that to
   * figure out rotational isotropicity
   */
  // TODO (first change the behavior of BondStereopermutators)
}

DirectedConformerGenerator::DirectedConformerGenerator(Molecule molecule)
  : _pImpl(std::make_unique<Impl>(std::move(molecule))) {}

DirectedConformerGenerator::DirectedConformerGenerator(DirectedConformerGenerator&& other) noexcept = default;
DirectedConformerGenerator& DirectedConformerGenerator::operator = (DirectedConformerGenerator&& other) noexcept = default;

DirectedConformerGenerator::~DirectedConformerGenerator() = default;

} // namespace molassembler
} // namespace Scine
