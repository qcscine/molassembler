#ifndef INCLUDE_TEMPLE_IDENTITY_H
#define INCLUDE_TEMPLE_IDENTITY_H

#include <utility>

/*!@file
 *
 * Provides the identity functor.
 */

namespace temple {

struct Identity {
  template<typename U>
  constexpr auto operator()(U&& v) const noexcept {
    return std::forward<U>(v);
  }
};

} // namespace temple

#endif
