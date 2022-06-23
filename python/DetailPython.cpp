/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "TypeCasters.h"
#include "pybind11/eigen.h"

#include "Molassembler/Detail/Cartesian.h"

void init_detail(pybind11::module& m) {
  using namespace Scine::Molassembler;

  auto detail = m.def_submodule("detail");
  detail.doc() = "Detail submodule";

  detail.def(
    "plane_of_best_fit_rmsd",
    &Cartesian::planeOfBestFitRmsd,
    pybind11::arg("positions"),
    pybind11::arg("indices"),
    "Fits a plane to indices and calculates its RMS deviation"
  );
}
