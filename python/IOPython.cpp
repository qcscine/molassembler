/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/eigen.h"

#include "Molassembler/IO.h"
#include "Molassembler/IO/SmilesParser.h"
#include "Molassembler/Molecule.h"

#include "Utils/Typenames.h"

void init_io(pybind11::module& m) {
  using namespace Scine::Molassembler;
  using namespace Scine::Utils;

  auto io = m.def_submodule("io");
  io.doc() = R"(IO Submodule)";

  auto experimental = io.def_submodule("experimental");
  experimental.doc() = R"(Experimental IO submodule)";

  experimental.def(
    "from_smiles_multiple",
    &IO::Experimental::parseSmiles,
    pybind11::arg("smiles_str"),
    R"delim(
      Parse a smiles string containing possibly multiple molecules

      The smiles parser is implemented according to the OpenSMILES spec. It
      supports the following features:
      - Arbitrarily many molecules in a string
      - Isotope markers (as long as they exist in Utils::ElementType)
      - Valence filling of the organic subset
      - Set shapes from VSEPR using (possibly) supplied charge
      - Ring closures
      - Stereo markers
        - Double bond
        - Tetrahedral (@ / @@ / @TH1 / @TH2)
        - Square planar (@SP1 - @SP3)
        - Trigonal bipyramidal (@TB1 - @TB20)
        - Octahedral (@OH1 - @OH30)

      :param smiles_str: A smiles string containing possibly multiple molecules
      :rtype: List[scine_molassembler.Molecule]

      >>> methane_and_ammonia = from_smiles_multiple("C.[NH4+]")
      >>> len(methane_and_ammonia) == 2
      True
    )delim"
  );

  experimental.def(
    "from_smiles",
    &IO::Experimental::parseSmilesSingleMolecule,
    pybind11::arg("smiles_str"),
    R"delim(
      Parse a smiles string containing only a single molecule

      The smiles parser is implemented according to the OpenSMILES spec. It
      supports the following features:
      - Arbitrarily many molecules in a string
      - Isotope markers (as long as they exist in Utils::ElementType)
      - Valence filling of the organic subset
      - Set shapes from VSEPR using (possibly) supplied charge
      - Ring closures
      - Stereo markers
        - Double bond
        - Tetrahedral (@ / @@ / @TH1 / @TH2)
        - Square planar (@SP1 - @SP3)
        - Trigonal bipyramidal (@TB1 - @TB20)
        - Octahedral (@OH1 - @OH30)

      :param smiles_str: A smiles string containing a single molecule
      :rtype: scine_molassembler.Molecule

      >>> import scine_utilities as utils
      >>> methane = from_smiles("C")
      >>> methane.graph.N == 4
      True
      >>> cobalt_complex = from_smiles("Br[Co@OH12](Cl)(I)(F)(S)C")
      >>> cobalt_index = cobalt_complex.graph.atoms_of_element(utils.ElementType.Co)[0]
      >>> permutator = cobalt_complex.stereopermutators.option(cobalt_index)
      >>> permutator is not None
      True
      >>> permutator.assigned is not None
      True
    )delim"
  );

  /* Line notations */
  pybind11::class_<IO::LineNotation> lineNotation(
    io,
    "LineNotation",
    "Generates :class:`Molecule` instances from line notations of molecules"
  );

  lineNotation.def_property_readonly_static(
    "enabled",
    [](pybind11::object /* self */) -> bool {
      return IO::LineNotation::enabled();
    },
    "Checks whether the ``obabel`` binary is found in your PATH"
  );

  lineNotation.def_static(
    "from_canonical_smiles",
    &IO::LineNotation::fromCanonicalSMILES,
    pybind11::arg("canonical_smiles"),
    "Construct a single :class:`Molecule` from a canonical SMILES string"
  );

  lineNotation.def_static(
    "from_isomeric_smiles",
    &IO::LineNotation::fromIsomericSMILES,
    pybind11::arg("isomeric_smiles"),
    "Construct a single :class:`Molecule` from an isomeric SMILES string"
  );

  lineNotation.def_static(
    "from_inchi",
    &IO::LineNotation::fromInChI,
    pybind11::arg("inchi"),
    "Construct a single :class:`Molecule` from an InChI string"
  );

  io.def(
    "read",
    &IO::read,
    pybind11::arg("filename"),
    R"delim(
      Reads a single :class:`Molecule` from a file. Interprets the file format from its
      extension. Supported formats:
      - mol: MOLFile V2000
      - xyz: XYZ file
      - cbor/bson/json: Serialization formats of molecules

      :param filename: File to read.
    )delim"
  );

  io.def(
    "split",
    &IO::split,
    pybind11::arg("filename"),
    R"delim(
      Reads multiple molecules from a file. Interprets the file format from its
      extension just like read(). Note that serializations of molecules contain
      only a single :class:`Molecule`. Use read() instead.

      :param filename: File to read.
    )delim"
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
    R"delim(
      Write a :class:`Molecule` and its positions to a file

      :param filename: File to write to. File format is interpreted from this
        parameter's file extension.
      :param molecule: :class:`Molecule` to write to file
      :param positions: Positions of molecule's atoms in bohr
    )delim"
  );

  io.def(
    "write",
    pybind11::overload_cast<const std::string&, const Molecule&>(&IO::write),
    pybind11::arg("filename"),
    pybind11::arg("molecule"),
    R"delim(
      Write a :class:`Molecule` serialization with the endings json/cbor/bson
      or a graph representation with ending dot/svg to a file.

      :param filename: File to write to. File format is interpreted from this
        parameter's file extension
      :param molecule: :class:`Molecule` to write to file
    )delim"
  );
}
