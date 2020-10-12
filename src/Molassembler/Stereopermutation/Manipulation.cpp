/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "Molassembler/Stereopermutation/Manipulation.h"

#include "Molassembler/Shapes/Data.h"
#include "Molassembler/Stereopermutation/RotationEnumerator.h"
#include "Molassembler/Temple/Functional.h"
#include "Molassembler/Temple/Permutations.h"

#include "boost/optional.hpp"
#include "boost/integer/common_factor_rt.hpp"

#include <unordered_map>

namespace Scine {
namespace Molassembler {
namespace Stereopermutations {

namespace {

inline unsigned gcd(const std::vector<unsigned>& c) {
  assert(!c.empty());
  auto iter = std::begin(c);
  const auto end = std::end(c);
  unsigned result = *iter;
  while(++iter != end) {
    result = boost::integer::gcd(result, *iter);
    if(result == 1) { return 1; }
  }

  return result;
}

} // namespace

inline void checkArguments(const Stereopermutation& s, const Shapes::Shape shape) {
  if(s.characters.size() != Shapes::size(shape)) {
    throw std::invalid_argument("Stereopermutation character count does not match shape size");
  }
}

std::vector<Stereopermutation> generateAllRotations(Stereopermutation s, const Shapes::Shape shape) {
  checkArguments(s, shape);
  RotationEnumerator enumerator {std::move(s), shape};
  return enumerator.all();
}

bool rotationallySuperimposable(Stereopermutation a, const Stereopermutation& b, const Shapes::Shape shape) {
  checkArguments(a, shape);
  checkArguments(b, shape);

  if(a == b) {
    return true;
  }

  RotationEnumerator enumerator {std::move(a), shape};

  while(auto rotationOption = enumerator.next()) {
    if(*rotationOption == b) {
      return true;
    }
  }

  return false;
}

boost::optional<bool> enantiomer(
  const Stereopermutation& a,
  const Stereopermutation& b,
  const Shapes::Shape shape
) {
  checkArguments(a, shape);
  checkArguments(b, shape);

  /* Generate the mirror image of *this and check whether it is rotationally
   * superimposable with other.
   */
  const auto& mirrorPermutation = Shapes::mirror(shape);

  /* If the mirror were to yield an identity permutation, it is represented
   * as an empty constexpr array (now a vector, so:)
   */
  if(mirrorPermutation.empty()) {
    return boost::none;
  }

  return rotationallySuperimposable(
    a.applyPermutation(mirrorPermutation),
    b,
    shape
  );
}

bool hasTransArrangedLinks(
  const Stereopermutation& s,
  const Shapes::Shape shape
) {
  checkArguments(s, shape);

  return Temple::any_of(
    s.links,
    [shape](const auto& link) {
      return Shapes::angleFunction(shape)(link.first, link.second) == M_PI;
    }
  );
}

Uniques uniques(
  const Stereopermutation& base,
  const Shapes::Shape shape,
  const bool removeTransSpanningGroups
) {
  checkArguments(base, shape);

  const unsigned S = Shapes::size(shape);
  auto permutation = Temple::iota<Shapes::Vertex>(S);
  auto stereopermutation = base;

  /* In case we want to skip trans pairs, the initial stereopermutation must also not
   * have any trans spanning pairs
   */
  if(removeTransSpanningGroups) {
    while(hasTransArrangedLinks(stereopermutation, shape)) {
      bool wasLastPermutation = !Temple::next_permutation(permutation);
      if(wasLastPermutation) {
        /* This can happen, e.g. in square-planar AAAB with
         * links: {0, 3}, {1, 3}, {2, 3}, every possible permutation contains
         * trans-arranged pairs. Then we return an empty vector.
         */
        return {};
      }

      stereopermutation = base.applyPermutation(permutation);
    }
  }

  auto minElement = [](const auto& c) {
    assert(!c.empty());
    return *std::min_element(std::begin(c), std::end(c));
  };

  // Generate all rotations of the base stereopermutation
  auto rotations = generateAllRotations(stereopermutation, shape);

  // The lowest rotation of the passed stereopermutation is the first unique stereopermutation
  Uniques unordered;
  unordered.list.push_back(minElement(rotations));
  unordered.weights.push_back(1);

  // Map of rotations to the index in the unordered so we can forward matches to those weight counters
  std::unordered_map<Stereopermutation, unsigned, boost::hash<Stereopermutation>> rotationCounterMap;
  for(const auto& rotation : rotations) {
    rotationCounterMap.emplace(rotation, 0);
  }

  // Go through all possible permutations of columns
  while(Temple::next_permutation(permutation)) {
    stereopermutation = base.applyPermutation(permutation);
    if(removeTransSpanningGroups && hasTransArrangedLinks(stereopermutation, shape)) {
      continue;
    }

    // Is the current stereopermutation not contained within the set of rotations?
    auto findIter = rotationCounterMap.find(stereopermutation);
    if(findIter == std::end(rotationCounterMap)) {
      // If so, it is a unique stereopermutation, generate all rotations
      rotations = generateAllRotations(stereopermutation, shape);
      unordered.list.push_back(minElement(rotations));
      unordered.weights.push_back(1);
      const unsigned key = unordered.weights.size() - 1;

      // Move from rotations
      for(auto& rotation : rotations) {
        rotationCounterMap.emplace(std::move(rotation), key);
      }
    } else {
      const unsigned weightIndex = findIter->second;
      ++unordered.weights.at(weightIndex);
    }
  }

  // Discover an ordering permutation
  const unsigned C = unordered.list.size();
  const std::vector<unsigned> order = Temple::sorted(
    Temple::iota<unsigned>(C),
    [&](const unsigned i, const unsigned j) -> bool {
      return unordered.list.at(i) < unordered.list.at(j);
    }
  );

  // Order the uniques using the discovered ordering permutation
  Uniques ordered;
  ordered.list.reserve(C);
  ordered.weights.resize(C);

  for(unsigned i = 0; i < C; ++i) {
    ordered.list.push_back(std::move(unordered.list.at(order.at(i))));
    ordered.weights.at(i) = unordered.weights.at(order.at(i));
  }

  // Divide the weights by their gcd
  const unsigned d = gcd(ordered.weights);
  for(unsigned& x : ordered.weights) {
    x /= d;
  }

  return ordered;
}

} // namespace Stereopermutations
} // namespace Molassembler
} // namespace Scine
