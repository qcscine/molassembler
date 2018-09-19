// Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
// See LICENSE.txt for details.

#ifndef INCLUDE_MOLASSEMBLER_TEMPLE_PREPROCESSOR_H
#define INCLUDE_MOLASSEMBLER_TEMPLE_PREPROCESSOR_H

/*! @file
 *
 * @brief Defines a set of useful preprocessor macros
 *
 * @note When applying the contained function attributes, remember that
 * exception branch correctness must be considered during e.g. common
 * subexpression elimination. It is helpful in many cases to additionally
 * annotate low-level often-used functions with the noexcept specifier (if
 * possible) to enable more optimizations.
 *
 * @note MSVC does not implement any of these attributes, so we must replace
 * them with empty defines if compiling with it.
 */

/*! Applies a weak purity function attribute
 *
 * The weak purity function attribute implies
 * - return value depends only on arguments and/or global variables (class
 *   member variables are implicitly an argument, as this* is captured)
 * - does not depend on volatile values
 * - does not cause side effects (and does not call any impure functions)
 * - is safe to call fewer times than the program says
 *
 * @warning This attribute is unchecked! Applying this attribute to functions
 * that do not adhere to the contract leads to undefined behavior.
 *
 * @note Do not expect compilers to aggressively take use of purity attributes
 * without also providing exception guarantees. See
 * https://kristerw.blogspot.com/2016/12/gcc-attributepure-and-c-exceptions.html
 */
#define PURITY_WEAK [[gnu::pure]]

/*! Applies a strong purity function attribute
 *
 * The strong purity function attribute implies all weak purity implications AND
 * - the return value depends only on function arguments
 * - its arguments are values, not pointers or references
 * - does not call non-pure or weak purity attributed functions
 *
 * @warning This attribute is unchecked! Applying this attribute to functions
 * that do not adhere to the contract leads to undefined behavior.
 *
 * @warning Be especially careful with classes. Simple getter functions that
 * return a member by value cannot be strongly pure as they refer to state that
 * is not encompassed by the function arguments. They may be weakly pure,
 * however.
 *
 * @warning Be REALLY careful with this attribute on template classes' members.
 * Make no unfounded assumptions about the template parameters behavior. If
 * needed, constrain the template parameter with SFINAE.
 *
 * @warning Any const class member function without arguments CANNOT be of
 * strong purity, since this is a pointer. Any access to class members within
 * the function violates the contract.
 *
 * @note Do not expect compilers to aggressively take use of purity attributes
 * without also providing exception guarantees. See
 * https://kristerw.blogspot.com/2016/12/gcc-attributepure-and-c-exceptions.html
 */
#define PURITY_STRONG [[gnu::const]]

#endif
