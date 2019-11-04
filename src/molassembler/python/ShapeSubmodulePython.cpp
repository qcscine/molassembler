/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "shapes/Data.h"

void init_shape_submodule(pybind11::module& m) {
  using namespace Scine;

  auto symmetrySubmodule = m.def_submodule("shape");
  symmetrySubmodule.doc() = R"(Shape submodule)";

  pybind11::enum_<Symmetry::Shape> shape(symmetrySubmodule, "Shape", "Shape enum");

  shape.value("Line", Symmetry::Shape::Line)
    .value("Bent", Symmetry::Shape::Bent)
    .value("EquilateralTriangle", Symmetry::Shape::EquilateralTriangle)
    .value("VacantTetrahedron", Symmetry::Shape::VacantTetrahedron)
    .value("T", Symmetry::Shape::T)
    .value("Tetrahedron", Symmetry::Shape::Tetrahedron)
    .value("Square", Symmetry::Shape::Square)
    .value("Seesaw", Symmetry::Shape::Seesaw)
    .value("TrigonalPyramid", Symmetry::Shape::TrigonalPyramid)
    .value("SquarePyramid", Symmetry::Shape::SquarePyramid)
    .value("TrigonalBipyramid", Symmetry::Shape::TrigonalBipyramid)
    .value("Pentagon", Symmetry::Shape::Pentagon)
    .value("Octahedron", Symmetry::Shape::Octahedron)
    .value("TrigonalPrism", Symmetry::Shape::TrigonalPrism)
    .value("PentagonalPyramid", Symmetry::Shape::PentagonalPyramid)
    .value("PentagonalBipyramid", Symmetry::Shape::PentagonalBipyramid)
    .value("SquareAntiprism", Symmetry::Shape::SquareAntiprism);

  shape.def("__str__", &Symmetry::name);

  symmetrySubmodule.def(
    "name_from_str",
    &Symmetry::nameFromString,
    pybind11::arg("name_str"),
    "Fetch a shape name from its string representation"
  );

  symmetrySubmodule.def(
    "size",
    &Symmetry::size,
    pybind11::arg("shape"),
    "Number of substituent positions in a shape"
  );
}
