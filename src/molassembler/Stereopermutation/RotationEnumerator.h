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
namespace Molassembler {
namespace Stereopermutations {

/*! @brief Enumerate rotations of stereopermutations in shapes
 *
 * Construction of class creates the first stereopermutation, repeated calls
 * to @p next generate successive rotations.
 *
 * Idea behind the algorithm and data structure is that starting from any
 * structure, we apply a rotation. If the result is new, we add it as a link in
 * a chain of structures. If the result has already been discovered, we try
 * again, but with a different rotation.
 *
 * Having written that down, it's probably not too difficult to rewrite this
 * into a recursive solution. This is more of a backtracking algorithm.
 */
class MASM_EXPORT RotationEnumerator {
public:
  using RotationSetType = std::vector<Stereopermutation>;

  //! Sets the initial state
  RotationEnumerator(Stereopermutation initial, Shapes::Shape s);

  /*! @brief Generates a new rotation or None if all rotations have been discovered
   *
   * For pattern `while(auto stereopermutationOption = enumerator.next())`
   */
  boost::optional<const Stereopermutation&> next();

  /*! @brief Yields set of all rotations
   *
   * This is just an accessor if all rotations have already been enumerated
   * with @p next, otherwise it first generates all rotations.
   */
  const RotationSetType& all();

private:
  bool incrementable() const;
  void increment();

  //! Link in chain
  struct Link {
    unsigned rotationIndex;
    Stereopermutation permutation;

    Link(unsigned i, Stereopermutation s);
  };

  Shapes::Shape shape;
  unsigned linkLimit;
  std::vector<Link> chain;
  RotationSetType rotations;
};

} // namespace Stereopermutations
} // namespace Molassembler
} // namespace Scine

#endif
