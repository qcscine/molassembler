/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "molassembler/Molecule.h"
#include "molassembler/Serialization.h"

/* Note: if you ever want to, you can get python objects from nlohmann's json
 * representation using the code from here:
 * https://github.com/pybind/pybind11/issues/1627
 *
 * Would require forward-declaring it in Serialization, though, and including
 * the interface library here (and linking it).
 */

std::vector<std::uint8_t> binaryFromPythonBytes(const pybind11::bytes& bytes) {
  std::string binaryStr = bytes;
  const char* binaryStrC = binaryStr.c_str();
  std::vector<std::uint8_t> binary;
  unsigned size = binaryStr.size();
  binary.resize(size);

  for(unsigned i = 0; i < size; ++i) {
    // ooh, pointer arithmetic, don't we just love that, clang-tidy?
    binary[i] = *(reinterpret_cast<const std::uint8_t*>(binaryStrC) + i);
  }
  return binary;
}

pybind11::bytes pythonBytesFromBinary(const std::vector<std::uint8_t>& binary) {
  return {
    reinterpret_cast<const char*>(&binary[0]),
    binary.size() * sizeof(std::uint8_t)
  };
}

pybind11::bytes JSONToBytes(
  const Scine::molassembler::JSONSerialization& serialization,
  const Scine::molassembler::JSONSerialization::BinaryFormat format
) {
  return pythonBytesFromBinary(serialization.toBinary(format));
}

Scine::molassembler::JSONSerialization JSONFromBytes(
  const pybind11::bytes& bytes,
  const Scine::molassembler::JSONSerialization::BinaryFormat format
) {
  return Scine::molassembler::JSONSerialization(binaryFromPythonBytes(bytes), format);
}

std::string toString(const Scine::molassembler::JSONSerialization& s) {
  return s;
}

void init_serialization(pybind11::module& m) {
  using namespace Scine::molassembler;

  pybind11::class_<JSONSerialization> serialization(
    m,
    "JSONSerialization",
    R"delim(
      Class representing a compact JSON serialization of a molecule

      >>> # Demonstrate a serialize-deserialize loop
      >>> import molassembler as masm
      >>> spiro = masm.io.experimental.from_smiles("C12(CCCC1)CCC2")
      >>> serializer = masm.JSONSerialization(spiro)
      >>> bson_format = masm.JSONSerialization.BinaryFormat.BSON
      >>> spiro_as_bson = serializer.to_binary(bson_format)
      >>> type(spiro_as_bson)
      <class 'bytes'>
      >>> bson_in_b64 = serializer.base_64_encode(spiro_as_bson)
      >>> type(bson_in_b64)
      <class 'str'>
      >>> reverted_bson = masm.JSONSerialization.base_64_decode(bson_in_b64)
      >>> serializer = masm.JSONSerialization(reverted_bson, bson_format)
      >>> reverted = serializer.to_molecule()
      >>> reverted == spiro # Compare the deserialized molecule
      True
    )delim"
  );

  pybind11::enum_<JSONSerialization::BinaryFormat> binaryFormat(
    serialization,
    "BinaryFormat",
    "Specifies the type of JSON binary format"
  );
  binaryFormat.value("CBOR", JSONSerialization::BinaryFormat::CBOR, "Compact Binary Object Representation")
   .value("BSON", JSONSerialization::BinaryFormat::BSON, "Binary JSON")
   .value("MsgPack", JSONSerialization::BinaryFormat::MsgPack, "MsgPack")
   .value("UBJSON", JSONSerialization::BinaryFormat::UBJSON, "Universal Binary JSON");

  serialization.def(
    pybind11::init<const std::string&>(),
    pybind11::arg("json_string"),
    "Parse a JSON molecule serialization"
  );

  serialization.def(
    pybind11::init<const Molecule&>(),
    pybind11::arg("molecule"),
    "Serialize a molecule into JSON"
  );

  serialization.def(
    pybind11::init(&JSONFromBytes),
    pybind11::arg("bytes"),
    pybind11::arg("binary_format"),
    "Parse a binary JSON molecule serialization"
  );

  serialization.def_static(
    "base_64_encode",
    [](const pybind11::bytes& bytes) -> std::string {
      return JSONSerialization::base64Encode(binaryFromPythonBytes(bytes));
    },
    pybind11::arg("binary"),
    "Encode binary format as base-64 string"
  );

  serialization.def_static(
    "base_64_decode",
    [](const std::string base64) -> pybind11::bytes {
      return pythonBytesFromBinary(JSONSerialization::base64Decode(base64));
    },
    pybind11::arg("base_64_string"),
    "Decode base-64 string into binary"
  );

  serialization.def(
    "__str__",
    &JSONSerialization::operator std::string,
    "Dump the JSON serialization into a string"
  );

  serialization.def(
    "to_string",
    &JSONSerialization::operator std::string,
    "Dump the JSON serialization into a string"
  );

  serialization.def(
    "standardize",
    &JSONSerialization::standardize,
    "Standardize the internal JSON serialization (only for canonical molecules)"
  );

  serialization.def(
    "to_molecule",
    &JSONSerialization::operator Scine::molassembler::Molecule,
    "Construct a molecule from the serialization"
  );

  serialization.def(
    "to_binary",
    &JSONToBytes,
    pybind11::arg("binary_format"),
    "Serialize a molecule into a binary format"
  );
}
