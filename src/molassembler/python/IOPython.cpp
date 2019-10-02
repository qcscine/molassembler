/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"

#include "molassembler/IO.h"
#include "molassembler/Molecule.h"

#include "Utils/Typenames.h"

void init_io(pybind11::module& m) {
  using namespace Scine::molassembler;
  using namespace Scine::Utils;

  auto io = m.def_submodule("io");
  io.doc() = R"(IO Submodule)";

  /* Line notations */
  pybind11::class_<IO::LineNotation> lineNotation(
    io,
    "LineNotation",
    "Generates Molecule instances from line notations of molecules"
  );

  lineNotation.def_property_readonly_static(
    "enabled",
    [](pybind11::object /* self */) -> bool {
      return IO::LineNotation::enabled();
    },
    "Checks whether the `obabel` binary is found in your PATH"
  );

  lineNotation.def_static(
    "from_canonical_smiles",
    &IO::LineNotation::fromCanonicalSMILES,
    pybind11::arg("canonical_smiles"),
    "Construct a single molecule from a canonical SMILES string"
  );

  lineNotation.def_static(
    "from_isomeric_smiles",
    &IO::LineNotation::fromIsomericSMILES,
    pybind11::arg("isomeric_smiles"),
    "Construct a single molecule from an isomeric SMILES string"
  );

  lineNotation.def_static(
    "from_inchi",
    &IO::LineNotation::fromInChI,
    pybind11::arg("inchi"),
    "Construct a single molecule from an InChI string"
  );

  io.def(
    "read",
    &IO::read,
    pybind11::arg("filename"),
    "Reads a single molecule from a file"
  );

  io.def(
    "split",
    &IO::split,
    pybind11::arg("filename"),
    "Reads multiple molecules from a file"
  );

  io.def(
    "write",
    pybind11::overload_cast<
      const std::string&,
      const Molecule&,
      const PositionCollection&
    >(&IO::write),
    pybind11::arg("filename"),
    pybind11::arg("molecule"),
    pybind11::arg("positions"),
    "Write a file from a molecule and positions"
  );

  io.def(
    "write",
    pybind11::overload_cast<const std::string&, const Molecule&>(&IO::write),
    pybind11::arg("filename"),
    pybind11::arg("molecule"),
    "Write a Molecule serialization with the endings json/cbor/bson to a file."
  );
}
