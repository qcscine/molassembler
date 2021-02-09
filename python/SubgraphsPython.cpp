/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */

#include "TypeCasters.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/Graph.h"
#include "Molassembler/Subgraphs.h"

using namespace Scine::Molassembler;

namespace {

template<typename T, typename BimapPybindType>
void init_bimap_side(BimapPybindType& bimap, const std::string& name) {
  pybind11::class_<T> side(bimap, name.c_str());

  using ValueType = typename T::value_type;
  pybind11::class_<ValueType> sideValue(side, "ValueType");

  sideValue.def("__repr__", [](const ValueType& v) -> std::string {
    return "(" + std::to_string(v.first) + ", " + std::to_string(v.second) + ")";
  });

  sideValue.def("__getitem__", [](const ValueType& v, const unsigned i) -> AtomIndex {
    if(i == 0) {
      return v.first;
    }

    if(i == 1) {
      return v.second;
    }

    throw std::out_of_range("Only two entries in ValueType");
  });

  sideValue.def("__len__", [](const ValueType& /* v */) { return 2; });

  side.def(
    "__iter__",
    [](const T& inst) {
      return pybind11::make_iterator(inst.begin(), inst.end());
    },
    pybind11::keep_alive<0, 1>()
  );
}

void init_bimap(pybind11::module& m) {
  pybind11::class_<Subgraphs::IndexMap> bimap(
    m,
    "Bimap",
    R"delim(
      Atom index bimap

      A bimap is a map type in which the stored relationship direction can
      easily be reversed. This particular bimap maps atom indices onto one
      another.

      The python bindings for the underlying boost::bimap are pretty limited.
      You can access the sides of the bimap and iterate through the mappings,
      and each mapping acts a little like a pair.

      You can get pure python types with e.g. ``[tuple(p) for p in bimap.left]``

      >>> neopentane = io.experimental.from_smiles("CC(C)(C)C")
      >>> methyl = io.experimental.from_smiles("[CH3]")
      >>> matches = complete(methyl, neopentane)
      >>> first_match = matches[0]
      >>> list(first_match.left)  # Mapping from methyl indices to neopentane
      [(0, 0), (1, 5), (2, 6), (3, 7)]
      >>> list(first_match.right)  # Mapping from neopentane to methyl indices
      [(0, 0), (5, 1), (6, 2), (7, 3)]
    )delim"
  );

  init_bimap_side<Subgraphs::IndexMap::left_map>(bimap, "Left");
  bimap.def_readonly("left", &Subgraphs::IndexMap::left, "Access stored relationships from the left");
  init_bimap_side<Subgraphs::IndexMap::right_map>(bimap, "Right");
  bimap.def_readonly("right", &Subgraphs::IndexMap::right, "Access stored relationships from the right");

  bimap.def("__len__", &Subgraphs::IndexMap::size);
}

} // namespace

void init_subgraphs(pybind11::module& m) {
  auto subgraphs = m.def_submodule("subgraphs");
  subgraphs.doc() = R"(Submodule for subgraph algorithms)";

  init_bimap(subgraphs);

  pybind11::enum_<Subgraphs::VertexStrictness> vertexStrictness(
    subgraphs,
    "VertexStrictness",
    "Matching strictness for vertices in subgraph matching"
  );

  vertexStrictness.value(
    "ElementType",
    Subgraphs::VertexStrictness::ElementType,
    "Element types of vertices must match"
  );

  pybind11::enum_<Subgraphs::EdgeStrictness> edgeStrictness(
    subgraphs,
    "EdgeStrictness",
    "Matching strictness for edges in subgraph matching"
  );

  edgeStrictness.value(
    "Topographic",
    Subgraphs::EdgeStrictness::Topographic,
    "No constraints are set upon edges besides topography"
  );

  edgeStrictness.value(
    "BondType",
    Subgraphs::EdgeStrictness::BondType,
    "Bond types between edges must match"
  );

  subgraphs.def(
    "complete",
    pybind11::overload_cast<
      const Graph&,
      const Graph&,
      Subgraphs::VertexStrictness,
      Subgraphs::EdgeStrictness
    >(&Subgraphs::complete),
    pybind11::arg("needle"),
    pybind11::arg("haystack"),
    pybind11::arg("vertex_strictness") = Subgraphs::VertexStrictness::ElementType,
    pybind11::arg("edge_strictness") = Subgraphs::EdgeStrictness::Topographic
  );

  subgraphs.def(
    "complete",
    pybind11::overload_cast<
      const Molecule&,
      const Molecule&,
      Subgraphs::VertexStrictness,
      Subgraphs::EdgeStrictness
    >(&Subgraphs::complete),
    pybind11::arg("needle"),
    pybind11::arg("haystack"),
    pybind11::arg("vertex_strictness") = Subgraphs::VertexStrictness::ElementType,
    pybind11::arg("edge_strictness") = Subgraphs::EdgeStrictness::Topographic
  );

  subgraphs.def(
    "maximum",
    pybind11::overload_cast<
      const Graph&,
      const Graph&,
      Subgraphs::VertexStrictness,
      Subgraphs::EdgeStrictness
    >(&Subgraphs::maximum),
    pybind11::arg("needle"),
    pybind11::arg("haystack"),
    pybind11::arg("vertex_strictness") = Subgraphs::VertexStrictness::ElementType,
    pybind11::arg("edge_strictness") = Subgraphs::EdgeStrictness::Topographic
  );

  subgraphs.def(
    "maximum",
    pybind11::overload_cast<
      const Molecule&,
      const Molecule&,
      Subgraphs::VertexStrictness,
      Subgraphs::EdgeStrictness
    >(&Subgraphs::maximum),
    pybind11::arg("needle"),
    pybind11::arg("haystack"),
    pybind11::arg("vertex_strictness") = Subgraphs::VertexStrictness::ElementType,
    pybind11::arg("edge_strictness") = Subgraphs::EdgeStrictness::Topographic
  );
}
