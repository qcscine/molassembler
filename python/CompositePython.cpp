/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Stereopermutation/Composites.h"

void init_composite(pybind11::module& m) {
  using namespace Scine::Molassembler;
  using namespace Stereopermutations;

  pybind11::class_<Composite> composite(
    m,
    "Composite",
    R"delim(
      Stereopermutation generating object for composites of two shapes
    )delim"
  );

  pybind11::class_<Composite::AngleGroup> angleGroup(
    composite,
    "AngleGroup",
    R"delim(
      A group of shape vertices at an angle from the fused position
    )delim"
  );

  angleGroup.def_readonly(
    "angle",
    &Composite::AngleGroup::angle,
    "Angle of the vertex group from the fused position"
  );

  angleGroup.def_readonly(
    "vertices",
    &Composite::AngleGroup::vertices,
    "Shape vertices comprising the group"
  );

  angleGroup.def_readonly(
    "isotropic",
    &Composite::AngleGroup::isotropic,
    "Whether the ranking characters indicate that this group of shape vertices is isotropic in this shape"
  );

  pybind11::class_<Composite::OrientationState> orientationState(
    composite,
    "OrientationState",
    R"delim(
      Orientation of a shape along a fused bond and reduced information on ranking
    )delim"
  );

  orientationState.def_readonly(
    "shape",
    &Composite::OrientationState::shape,
    "Shape at this side of the bond"
  );

  orientationState.def_readonly(
    "fused_vertex",
    &Composite::OrientationState::fusedVertex,
    "Shape vertex that is fused to the other shape"
  );

  orientationState.def_readonly(
    "characters",
    &Composite::OrientationState::characters,
    "Ranking characters of all vertices of the shape"
  );

  orientationState.def_readonly(
    "identifier",
    &Composite::OrientationState::identifier,
    "An identifier to the shape source"
  );

  orientationState.def_property_readonly(
    "smallest_angle_group",
    &Composite::OrientationState::smallestAngleGroup,
    "Collects all coplanar indices that are closest to the fused shape vertex"
  );

  composite.def_property_readonly(
    "isotropic",
    &Composite::isIsotropic,
    "Whether the Composite is isotropic overall"
  );

  composite.def_property_readonly(
    "order",
    &Composite::order,
    "The higher number of relevant vertices of both sides"
  );

  composite.def_property_readonly(
    "orientations",
    [](const Composite& c) {
      return std::make_pair(
        c.orientations().first,
        c.orientations().second
      );
    },
    "Orientations of the composite"
  );

  composite.def(
    "__iter__",
    [](const Composite& c) { return pybind11::make_iterator(c.begin(), c.end()); },
    "Iterate through stereopermutations"
  );

  composite.def(
    "__len__",
    [](const Composite& c) { return c.permutations(); }
  );

  composite.def(
    "__getitem__",
    [](const Composite& c, const unsigned i) {
      return c.dihedrals(i);
    }
  );
}
