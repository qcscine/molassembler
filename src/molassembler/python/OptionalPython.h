/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#ifndef INCLUDE_MOLASSEMBLER_PYTHON_OPTIONAL_H
#define INCLUDE_MOLASSEMBLER_PYTHON_OPTIONAL_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "boost/optional.hpp"

namespace pybind11 {
namespace detail {

template <typename T>
struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};

} // namespace detail
} // namespace pybind11

#endif
