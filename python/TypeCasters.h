/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#ifndef INCLUDE_MOLASSEMBLER_PYTHON_TYPE_CASTERS_H
#define INCLUDE_MOLASSEMBLER_PYTHON_TYPE_CASTERS_H

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "boost/variant.hpp"
#include "boost/optional.hpp"
#include "molassembler/RankingInformation.h"

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

template <typename T>
struct type_caster<boost::optional<T>> : optional_caster<boost::optional<T>> {};

/* Type caster for SiteIndex <-> number */
template<>
struct type_caster<Scine::Molassembler::SiteIndex> {
  PYBIND11_TYPE_CASTER(Scine::Molassembler::SiteIndex, _("SiteIndex"));

  bool load(handle src, bool) {
    PyObject* source = src.ptr();
    // Try to convert the handle to an integer
    PyObject* tmp = PyNumber_Index(source);
    if(!tmp) {
      return false;
    }

    PyObject* exc = nullptr;
    const auto pyIndex = PyNumber_AsSsize_t(tmp, exc);
    value = Scine::Molassembler::SiteIndex(pyIndex);
    Py_DECREF(tmp);
    if(exc) {
      return false;
    }

    return true;
  }

  static handle cast(
    Scine::Molassembler::SiteIndex src,
    return_value_policy /* policy */,
    handle /* parent */
  ) {
    return PyLong_FromLong(src);
  }
};

} // namespace detail
} // namespace pybind11

#endif
