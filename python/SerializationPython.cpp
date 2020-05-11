/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"

#include "Molassembler/Molecule.h"
#include "Molassembler/Serialization.h"

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

pybind11::bytes JsonToBytes(
  const Scine::Molassembler::JsonSerialization& serialization,
  const Scine::Molassembler::JsonSerialization::BinaryFormat format
) {
  return pythonBytesFromBinary(serialization.toBinary(format));
}

Scine::Molassembler::JsonSerialization JsonFromBytes(
  const pybind11::bytes& bytes,
  const Scine::Molassembler::JsonSerialization::BinaryFormat format
) {
  return Scine::Molassembler::JsonSerialization(binaryFromPythonBytes(bytes), format);
}

std::string toString(const Scine::Molassembler::JsonSerialization& s) {
  return s;
}

void init_serialization(pybind11::module& m) {
  using namespace Scine::Molassembler;

  pybind11::class_<JsonSerialization> serialization(
    m,
    "JsonSerialization",
    R"delim(
      Class representing a compact JSON serialization of a molecule

      >>> # Demonstrate a serialize-deserialize loop
      >>> spiro = io.experimental.from_smiles("C12(CCCC1)CCC2")
      >>> serializer = JsonSerialization(spiro)
      >>> bson_format = JsonSerialization.BinaryFormat.BSON
      >>> spiro_as_bson = serializer.to_binary(bson_format)
      >>> bson_in_b64 = serializer.base_64_encode(spiro_as_bson)
      >>> reverted_bson = JsonSerialization.base_64_decode(bson_in_b64)
      >>> serializer = JsonSerialization(reverted_bson, bson_format)
      >>> reverted = serializer.to_molecule()
      >>> reverted == spiro # Compare the deserialized molecule
      True
    )delim"
  );

  pybind11::enum_<JsonSerialization::BinaryFormat> binaryFormat(
    serialization,
    "BinaryFormat",
    "Specifies the type of JSON binary format"
  );
  binaryFormat.value("CBOR", JsonSerialization::BinaryFormat::CBOR, "Compact Binary Object Representation")
   .value("BSON", JsonSerialization::BinaryFormat::BSON, "Binary JSON")
   .value("MsgPack", JsonSerialization::BinaryFormat::MsgPack, "MsgPack")
   .value("UBJSON", JsonSerialization::BinaryFormat::UBJSON, "Universal Binary JSON");

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
    pybind11::init(&JsonFromBytes),
    pybind11::arg("bytes"),
    pybind11::arg("binary_format"),
    "Parse a binary JSON molecule serialization"
  );

  serialization.def_static(
    "base_64_encode",
    [](const pybind11::bytes& bytes) -> std::string {
      return JsonSerialization::base64Encode(binaryFromPythonBytes(bytes));
    },
    pybind11::arg("binary"),
    "Encode binary format as base-64 string"
  );

  serialization.def_static(
    "base_64_decode",
    [](const std::string base64) -> pybind11::bytes {
      return pythonBytesFromBinary(JsonSerialization::base64Decode(base64));
    },
    pybind11::arg("base_64_string"),
    "Decode base-64 string into binary"
  );

  serialization.def(
    "__str__",
    &JsonSerialization::operator std::string,
    "Dump the JSON serialization into a string"
  );

  serialization.def(
    "to_string",
    &JsonSerialization::operator std::string,
    "Dump the JSON serialization into a string"
  );

  serialization.def(
    "standardize",
    &JsonSerialization::standardize,
    "Standardize the internal JSON serialization (only for canonical molecules)"
  );

  serialization.def(
    "to_molecule",
    &JsonSerialization::operator Scine::Molassembler::Molecule,
    "Construct a molecule from the serialization"
  );

  serialization.def(
    "to_binary",
    &JsonToBytes,
    pybind11::arg("binary_format"),
    "Serialize a molecule into a binary format"
  );
}
