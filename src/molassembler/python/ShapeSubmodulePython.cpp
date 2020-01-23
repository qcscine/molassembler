/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "shapes/Data.h"

void init_shape_submodule(pybind11::module& m) {
  using namespace Scine;

  auto symmetrySubmodule = m.def_submodule("shapes");
  symmetrySubmodule.doc() = R"(Shape submodule)";

  pybind11::enum_<shapes::Shape> shape(symmetrySubmodule, "Shape", "Shape enum");

  shape.value("Line", shapes::Shape::Line)
    .value("Bent", shapes::Shape::Bent)
    .value("EquilateralTriangle", shapes::Shape::EquilateralTriangle)
    .value("VacantTetrahedron", shapes::Shape::VacantTetrahedron)
    .value("T", shapes::Shape::T)
    .value("Tetrahedron", shapes::Shape::Tetrahedron)
    .value("Square", shapes::Shape::Square)
    .value("Seesaw", shapes::Shape::Seesaw)
    .value("TrigonalPyramid", shapes::Shape::TrigonalPyramid)
    .value("SquarePyramid", shapes::Shape::SquarePyramid)
    .value("TrigonalBipyramid", shapes::Shape::TrigonalBipyramid)
    .value("Pentagon", shapes::Shape::Pentagon)
    .value("Octahedron", shapes::Shape::Octahedron)
    .value("TrigonalPrism", shapes::Shape::TrigonalPrism)
    .value("PentagonalPyramid", shapes::Shape::PentagonalPyramid)
    .value("Hexagon", shapes::Shape::Hexagon)
    .value("PentagonalBipyramid", shapes::Shape::PentagonalBipyramid)
    .value("CappedOctahedron", shapes::Shape::CappedOctahedron)
    .value("CappedTrigonalPrism", shapes::Shape::CappedTrigonalPrism)
    .value("SquareAntiprism", shapes::Shape::SquareAntiprism)
    .value("Cube", shapes::Shape::Cube)
    .value("TrigonalDodecahedron", shapes::Shape::TrigonalDodecahedron)
    .value("HexagonalBipyramid", shapes::Shape::HexagonalBipyramid)
    .value("TricappedTrigonalPrism", shapes::Shape::TricappedTrigonalPrism)
    .value("CappedSquareAntiprism", shapes::Shape::CappedSquareAntiprism)
    .value("HeptagonalBipyramid", shapes::Shape::HeptagonalBipyramid)
    .value("BicappedSquareAntiprism", shapes::Shape::BicappedSquareAntiprism)
    .value("EdgeContractedIcosahedron", shapes::Shape::EdgeContractedIcosahedron)
    .value("Icosahedron", shapes::Shape::Icosahedron)
    .value("Cuboctahedron", shapes::Shape::Cuboctahedron);

  shape.def("__str__", &shapes::name);

  symmetrySubmodule.def(
    "name_from_str",
    &shapes::nameFromString,
    pybind11::arg("name_str"),
    "Fetch a shape name from its string representation"
  );

  symmetrySubmodule.def(
    "size",
    &shapes::size,
    pybind11::arg("shape"),
    "Number of substituent positions in a shape"
  );
}
