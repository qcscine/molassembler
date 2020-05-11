/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "molassembler/Shapes/Data.h"

void init_shape_submodule(pybind11::module& m) {
  using namespace Scine;

  auto symmetrySubmodule = m.def_submodule("shapes");
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
    .value("Hexagon", Shapes::Shape::Hexagon)
    .value("PentagonalBipyramid", Shapes::Shape::PentagonalBipyramid)
    .value("CappedOctahedron", Shapes::Shape::CappedOctahedron)
    .value("CappedTrigonalPrism", Shapes::Shape::CappedTrigonalPrism)
    .value("SquareAntiprism", Shapes::Shape::SquareAntiprism)
    .value("Cube", Shapes::Shape::Cube)
    .value("TrigonalDodecahedron", Shapes::Shape::TrigonalDodecahedron)
    .value("HexagonalBipyramid", Shapes::Shape::HexagonalBipyramid)
    .value("TricappedTrigonalPrism", Shapes::Shape::TricappedTrigonalPrism)
    .value("CappedSquareAntiprism", Shapes::Shape::CappedSquareAntiprism)
    .value("HeptagonalBipyramid", Shapes::Shape::HeptagonalBipyramid)
    .value("BicappedSquareAntiprism", Shapes::Shape::BicappedSquareAntiprism)
    .value("EdgeContractedIcosahedron", Shapes::Shape::EdgeContractedIcosahedron)
    .value("Icosahedron", Shapes::Shape::Icosahedron)
    .value("Cuboctahedron", Shapes::Shape::Cuboctahedron);

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
