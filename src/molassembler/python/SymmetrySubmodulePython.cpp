/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "chemical_symmetries/Symmetries.h"

void init_symmetry_submodule(pybind11::module& m) {
  using namespace Scine;

  auto symmetrySubmodule = m.def_submodule("symmetry");

  pybind11::enum_<Symmetry::Name> name(symmetrySubmodule, "SymmetryName", "Symmetry name enum");

  name.value("Linear", Symmetry::Name::Linear)
    .value("Bent", Symmetry::Name::Bent)
    .value("TrigonalPlanar", Symmetry::Name::TrigonalPlanar)
    .value("CutTetrahedral", Symmetry::Name::CutTetrahedral)
    .value("TShaped", Symmetry::Name::TShaped)
    .value("Tetrahedral", Symmetry::Name::Tetrahedral)
    .value("SquarePlanar", Symmetry::Name::SquarePlanar)
    .value("Seesaw", Symmetry::Name::Seesaw)
    .value("TrigonalPyramidal", Symmetry::Name::TrigonalPyramidal)
    .value("SquarePyramidal", Symmetry::Name::SquarePyramidal)
    .value("TrigonalBiPyramidal", Symmetry::Name::TrigonalBiPyramidal)
    .value("PentagonalPlanar", Symmetry::Name::PentagonalPlanar)
    .value("Octahedral", Symmetry::Name::Octahedral)
    .value("TrigonalPrismatic", Symmetry::Name::TrigonalPrismatic)
    .value("PentagonalPyramidal", Symmetry::Name::PentagonalPyramidal)
    .value("PentagonalBiPyramidal", Symmetry::Name::PentagonalBiPyramidal)
    .value("SquareAntiPrismatic", Symmetry::Name::SquareAntiPrismatic);

  name.def("__str__", &Symmetry::name);

  symmetrySubmodule.def("name_from_str", &Symmetry::nameFromString, "Fetch a symmetry name from its string representation");

  symmetrySubmodule.def("size", &Symmetry::size, "Number of substituent positions in a symmetry");
}
