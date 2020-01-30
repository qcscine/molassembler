/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
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
struct type_caster<Scine::molassembler::SiteIndex> {
  PYBIND11_TYPE_CASTER(Scine::molassembler::SiteIndex, _("SiteIndex"));

  bool load(handle src, bool) {
    PyObject* source = src.ptr();
    PyObject* tmp = PyNumber_Long(source);
    if(!tmp) {
      return false;
    }

    const long signedValue = PyLong_AsLong(tmp);
    value = Scine::molassembler::SiteIndex(signedValue);
    Py_DECREF(tmp);
    return !(signedValue == -1 && !PyErr_Occurred());
  }

  static handle cast(
    Scine::molassembler::SiteIndex src,
    return_value_policy /* policy */,
    handle /* parent */
  ) {
    return PyLong_FromLong(src);
  }
};

} // namespace detail
} // namespace pybind11

#endif
