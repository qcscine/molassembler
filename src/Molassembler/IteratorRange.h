/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Range type
 */

#ifndef INCLUDE_MOLASSEMBLER_RANGE_H
#define INCLUDE_MOLASSEMBLER_RANGE_H

#include "Molassembler/Export.h"
#include <utility>

namespace Scine {
namespace Molassembler {

/**
 * @brief Homogeneous pair of iterators with begin and end member fns
 * @tparam Iter Iterator type to store
 *
 * A function returning two iterators as a semantic range can return this to
 * make the function range-for compatible. Semantically, this class is little
 * different from a std::pair with homogeneous template arguments.
 */
template<typename Iter>
struct MASM_EXPORT IteratorRange {
  Iter first;
  Iter second;

  IteratorRange() = default;
  IteratorRange(Iter a, Iter b) : first(std::move(a)), second(std::move(b)) {}

  inline Iter begin() const {
    return first;
  }

  inline Iter end() const {
    return second;
  }
};

} // namespace Molassembler
} // namespace Scine

#endif
