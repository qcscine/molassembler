/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */

#include "stereopermutation/RotationEnumerator.h"

#include "shapes/Data.h"

namespace Scine {
namespace stereopermutation {

RotationEnumerator::Link::Link(unsigned i, Stereopermutation s) : rotationIndex(i), permutation(std::move(s)) {}

RotationEnumerator::RotationEnumerator(Stereopermutation initial, const shapes::Shape s)
  : shape(s),
    linkLimit(shapes::rotations(shape).size())
{
  chain.emplace_back(
    0u,
    initial
  );

  rotations.push_back(std::move(initial));
}

inline bool contains(const std::vector<Stereopermutation>& c, const Stereopermutation& s) {
  return std::find(std::begin(c), std::end(c), s) != std::end(c);
}

bool RotationEnumerator::incrementable() const {
  return chain.front().rotationIndex < linkLimit;
}

void RotationEnumerator::increment() {
  unsigned& lastChainRotationIndex = chain.back().rotationIndex;
  if(lastChainRotationIndex < linkLimit - 1) {
    ++lastChainRotationIndex;
  } else {
    while(chain.size() > 1 && chain.back().rotationIndex == linkLimit - 1) {
      chain.pop_back();
    }

    ++chain.back().rotationIndex;
  }
}

boost::optional<const Stereopermutation&> RotationEnumerator::next() {
  while(incrementable()) {
    const Link& lastLink = chain.back();
    Stereopermutation rotation = lastLink.permutation.applyPermutation(
      shapes::rotations(shape).at(lastLink.rotationIndex)
    );

    if(!contains(rotations, rotation)) {
      rotations.push_back(rotation);
      chain.emplace_back(
        0,
        std::move(rotation)
      );

      return chain.back().permutation;
    } else {
      increment();
    }
  }

  return boost::none;
}

const RotationEnumerator::RotationSetType& RotationEnumerator::all() {
  while(incrementable()) {
    const Link& lastLink = chain.back();
    Stereopermutation rotation = lastLink.permutation.applyPermutation(
      shapes::rotations(shape).at(lastLink.rotationIndex)
    );

    if(!contains(rotations, rotation)) {
      rotations.push_back(rotation);
      chain.emplace_back(
        0,
        std::move(rotation)
      );
    } else {
      increment();
    }
  }

  return rotations;
}

} // namespace stereopermutation
} // namespace Scine
