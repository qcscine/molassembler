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

  pybind11::enum_<Shapes::Shape> shape(symmetrySubmodule, "Shape", "Shape enum");

  shape.value("Line", Shapes::Shape::Line)
    .value("Bent", Shapes::Shape::Bent)
    .value("EquilateralTriangle", Shapes::Shape::EquilateralTriangle)
    .value("VacantTetrahedron", Shapes::Shape::VacantTetrahedron)
    .value("T", Shapes::Shape::T)
    .value("Tetrahedron", Shapes::Shape::Tetrahedron)
    .value("Square", Shapes::Shape::Square)
    .value("Seesaw", Shapes::Shape::Seesaw)
    .value("TrigonalPyramid", Shapes::Shape::TrigonalPyramid)
    .value("SquarePyramid", Shapes::Shape::SquarePyramid)
    .value("TrigonalBipyramid", Shapes::Shape::TrigonalBipyramid)
    .value("Pentagon", Shapes::Shape::Pentagon)
    .value("Octahedron", Shapes::Shape::Octahedron)
    .value("TrigonalPrism", Shapes::Shape::TrigonalPrism)
    .value("PentagonalPyramid", Shapes::Shape::PentagonalPyramid)
    .value("PentagonalBipyramid", Shapes::Shape::PentagonalBipyramid)
    .value("SquareAntiprism", Shapes::Shape::SquareAntiprism);

  shape.def("__str__", &Shapes::name);

  symmetrySubmodule.def(
    "name_from_str",
    &Shapes::nameFromString,
    pybind11::arg("name_str"),
    "Fetch a shape name from its string representation"
  );

  symmetrySubmodule.def(
    "size",
    &Shapes::size,
    pybind11::arg("shape"),
    "Number of substituent positions in a shape"
  );
}
