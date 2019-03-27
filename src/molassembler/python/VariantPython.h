/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#ifndef INCLUDE_MOLASSEMBLER_PYTHON_VARIANT_H
#define INCLUDE_MOLASSEMBLER_PYTHON_VARIANT_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "boost/variant.hpp"

namespace pybind11 {
namespace detail {

template <typename... Ts>
struct type_caster<boost::variant<Ts...>> : variant_caster<boost::variant<Ts...>> {};

template <>
struct visit_helper<boost::variant> {
    template <typename... Args>
    static auto call(Args &&...args) -> decltype(boost::apply_visitor(args...)) {
        return boost::apply_visitor(args...);
    }
};

} // namespace detail
} // namespace pybind11

#endif
