/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "pybind11/eigen.h"

#include "Molassembler/IO.h"
#include "Molassembler/IO/SmilesParser.h"
#include "Molassembler/IO/SmilesEmitter.h"
#include "Molassembler/Molecule.h"

#include "Utils/Typenames.h"

void init_io(pybind11::module& m) {
  using namespace Scine::Molassembler;
  using namespace Scine::Utils;

  auto io = m.def_submodule("io");
  io.doc() = R"(IO Submodule)";

  auto experimental = io.def_submodule("experimental");
  experimental.doc() = R"(
    Experimental
    ------------

    :note: Functions in this module are unstable and should be used with
      caution. Check your results. Upon stabilization, functions will be
      deprecated and move to a different module.
  )";

  experimental.def(
    "from_smiles_multiple",
    &IO::Experimental::parseSmiles,
    pybind11::arg("smiles_str"),
    R"delim(
      Parse a smiles string containing possibly multiple molecules

      The smiles parser is implemented according to the OpenSMILES spec. It
      supports the following features:

      - Arbitrarily many molecules in a string
      - Isotope markers
      - Valence filling of the organic subset
      - Set shapes from VSEPR
      - Ring closures
      - Stereo markers
        - Double bond
        - Tetrahedral
        - Square planar
        - Trigonal bipyramidal
        - Octahedral

      :param smiles_str: A smiles string containing possibly multiple molecules

      >>> methane_and_ammonia = io.experimental.from_smiles_multiple("C.[NH4+]")
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

      - Isotope markers
      - Valence filling of the organic subset
      - Set local shapes from VSEPR
      - Ring closures (and concatenation between dot-separated components)
      - Stereo markers
        - Double bond
        - Tetrahedral
        - Square planar
        - Trigonal bipyramidal
        - Octahedral

      :param smiles_str: A smiles string containing a single molecule

      >>> import scine_utilities as utils
      >>> methane = io.experimental.from_smiles("C")
      >>> methane.graph.V
      5
      >>> cobalt_complex = io.experimental.from_smiles("Br[Co@OH12](Cl)(I)(F)(S)C")
      >>> cobalt_index = cobalt_complex.graph.atoms_of_element(utils.ElementType.Co)[0]
      >>> permutator = cobalt_complex.stereopermutators.option(cobalt_index)
      >>> permutator is not None
      True
      >>> permutator.assigned is not None
      True
    )delim"
  );

  experimental.def(
    "emit_smiles",
    &IO::Experimental::emitSmiles,
    pybind11::arg("molecule"),
    R"delim(
      Generate a smiles string for a molecule

      :param molecule: Molecule to generate smiles string for
      :returns: A (partially) normalized openSMILES-standard compliant
        smiles string.

      :warning: This is a lossy serialization format! The openSMILES
        standard does not contain stereodescriptors for shapes other than
        the tetrahedron, square, trigonal bipyramid and octahedron.
        Generated smiles containing stereocenters with other shapes will
        not contain stereodescriptors for these centers.

      :note: Missing normalization: Aromaticity detection in kekulized
        graph to aromatic atom types.

      >>> biphenyl = io.experimental.from_smiles("c1ccccc1-c2ccccc2")
      >>> io.experimental.emit_smiles(biphenyl)
      'c1ccccc1-c2ccccc2'
    )delim"
  );

  /* Line notations */
  pybind11::class_<IO::LineNotation> lineNotation(
    io,
    "LineNotation",
    R"delim(
      Generates :class:`Molecule` instances from line notations of molecules
      via OpenBabel, if found in the runtime path.
    )delim"
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
