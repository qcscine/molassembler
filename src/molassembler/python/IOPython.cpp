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
    "Construct a single molecule from a canonical SMILES string"
  );

  lineNotation.def_static(
    "from_isomeric_smiles",
    &IO::LineNotation::fromIsomericSMILES,
    "Construct a single molecule from an isomeric SMILES string"
  );

  lineNotation.def_static(
    "from_inchi",
    &IO::LineNotation::fromInChI,
    "Construct a single molecule from an InChI string"
  );

  io.def(
    "read",
    &IO::read,
    "Reads a single molecule from a file"
  );

  io.def(
    "split",
    &IO::split,
    "Reads multiple molecules from a file"
  );

  io.def(
    "write",
    pybind11::overload_cast<
      const std::string&,
      const Molecule&,
      const PositionCollection&
    >(&IO::write),
    "Write a file from a molecule and positions"
  );

  io.def(
    "write",
    pybind11::overload_cast<const std::string&, const Molecule&>(&IO::write),
    "Write a Molecule serialization (JSON / CBOR) directly to a file."
  );
}
