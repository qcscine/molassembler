/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Helper to enumerate rotations of stereopermutations in shapes
 */

#ifndef INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_ROTATION_ENUMERATOR_H
#define INCLUDE_MOLASSEMBLER_STEREOPERMUTATIONS_ROTATION_ENUMERATOR_H

#include "molassembler/Stereopermutation/Stereopermutation.h"
#include "molassembler/Shapes/Shapes.h"
#include "boost/functional/hash.hpp"
#include "boost/optional/optional_fwd.hpp"

#include <unordered_set>

namespace Scine {
namespace stereopermutation {

class RotationEnumerator {
public:
  using RotationSetType = std::vector<Stereopermutation>;

  RotationEnumerator(Stereopermutation initial, shapes::Shape s);

  boost::optional<const Stereopermutation&> next();
  const RotationSetType& all();

private:
  bool incrementable() const;
  void increment();

  struct Link {
    unsigned rotationIndex;
    Stereopermutation permutation;

    Link(unsigned i, Stereopermutation s);
  };

  shapes::Shape shape;
  unsigned linkLimit;
  std::vector<Link> chain;
  RotationSetType rotations;
};

} // namespace stereopermutation
} // namespace Scine

#endif
