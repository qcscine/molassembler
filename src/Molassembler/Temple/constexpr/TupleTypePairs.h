/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 * @brief Pass pairs of types contained in a tuple type to template functions
 */

#ifndef INCLUDE_TEMPLE_TUPLE_TYPE_PAIRS_H
#define INCLUDE_TEMPLE_TUPLE_TYPE_PAIRS_H

#include "Molassembler/Temple/constexpr/UpperTriangularMatrix.h"
#include "Molassembler/Temple/constexpr/TupleType.h"

namespace Scine {
namespace Molassembler {
namespace Temple {
namespace Tuples {
namespace Detail {

/*!
 * Implementation of a special case of mapping in which all possible pairs of
 * types contained in the supplied tuple are passed to a template function and
 * evaluated.
 */
template<
  typename TupleType,
  template<typename ...> class TemplateFunction,
  size_t ... Inds
> constexpr auto mapAllPairsHelper(
  std::index_sequence<Inds...> /* inds */
) {
  constexpr size_t N = std::tuple_size<TupleType>::value;

  /* Although you might think this is awful, not immediately consuming the
   * generated index pairs, there is basically no elegant way to generate and
   * consume in two places without incurring another function that needs to be
   * instantiated for every pair of indices separately, which would be worse.
   */
  constexpr auto indexPairs = makeArray(
    UpperTriangularMatrixImpl::index_conversion::toDoubleIndex<N>(Inds)...
  );

  return makeArray(
    handleValueVariants<
      TemplateFunction<
        std::tuple_element_t<indexPairs.at(Inds).first, TupleType>,
        std::tuple_element_t<indexPairs.at(Inds).second, TupleType>
      >
    >()...
  );
}

} // namespace Detail

/*! @brief Map all pairs of types in a tuple typedefs into a template function
 *
 * Takes a tuple type and a template function that accepts pairs of the tuple's
 * contained types. Returns an array containing the result of the mapping for
 * all distinct possible type pairs.
 *
 * With distinct, we mean that pairs of the same types are not created unless a
 * type is explicitly repeated in the tuple type: Say the tuple contains Apple,
 * Banana and Orange, then the generated pairs are: Apple & Banana,
 * Apple & Orange and Banana & Orange.
 *
 * @complexity{@math{\Theta(N^2)}}
 *
 * @tparam TupleType tuple typedef containing the types whose pairs should be
 *   mapped
 * @tparam TemplateFunction template metafunction with a static constexpr value
 *   member or member function
 */
template<
  typename TupleType,
  template<typename ...> class TemplateFunction
> constexpr auto mapAllPairs() {
  constexpr size_t N = std::tuple_size<TupleType>::value;

  return Detail::mapAllPairsHelper<TupleType, TemplateFunction>(
    std::make_index_sequence<(N * N - N) / 2>()
  );
}

} // namespace Tuples
} // namespace Temple
} // namespace Molassembler
} // namespace Scine

#endif
