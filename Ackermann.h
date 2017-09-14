#ifndef INCLUDE_CONSTEXPR_MAGIC_ACKERMANN_H
#define INCLUDE_CONSTEXPR_MAGIC_ACKERMANN_H

/*! @file
 *
 * Defines the ackermann function if you ever want to incur heavy computational
 * cost in any particular place.
 */

namespace ConstexprMagic {

/*!
 * The Ackermann-Peter function, whose values increase very rapidly, and where
 * the computation of (4, 1) can take quite a bit of time and memory.
 */
constexpr unsigned ackermann(const unsigned& m, const unsigned& n) {
  if (m == 0) {
    return n + 1;
  }

  if (n == 0) {
    return ackermann(m - 1, 1);
  }

  return ackermann(
    m - 1,
    ackermann(m, n - 1)
  );
}

} // namespace ConstexprMagic

#endif
