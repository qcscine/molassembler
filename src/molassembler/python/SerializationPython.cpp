/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "molassembler/Molecule.h"
#include "molassembler/Serialization.h"

std::string toJSON(const Scine::molassembler::Molecule& mol) {
  return Scine::molassembler::JSONSerialization(mol);
}

Scine::molassembler::Molecule fromJSON(const std::string& jsonStr) {
  return Scine::molassembler::JSONSerialization(jsonStr);
}

pybind11::bytes toCBORBytes(const Scine::molassembler::Molecule& mol) {
  auto binary = Scine::molassembler::JSONSerialization(mol).toBinary(
    Scine::molassembler::JSONSerialization::BinaryFormat::CBOR
  );

  return {
    reinterpret_cast<const char*>(&binary[0]),
    binary.size() * sizeof(std::uint8_t)
  };
}

pybind11::bytes toBSONBytes(const Scine::molassembler::Molecule& mol) {
  auto binary = Scine::molassembler::JSONSerialization(mol).toBinary(
    Scine::molassembler::JSONSerialization::BinaryFormat::BSON
  );

  return {
    reinterpret_cast<const char*>(&binary[0]),
    binary.size() * sizeof(std::uint8_t)
  };
}

std::vector<std::uint8_t> binaryFromPybind11Bytes(const pybind11::bytes& bytes) {
  std::string binaryStr = bytes;
  const char* binaryStrC = binaryStr.c_str();

  std::vector<std::uint8_t> binary;
  unsigned size = binaryStr.size();
  binary.resize(size);
  for(unsigned i = 0; i < size; ++i) {
    binary[i] = *reinterpret_cast<std::uint8_t*>(binaryStrC[i]);
  }

  return binary;
}

Scine::molassembler::Molecule fromCBORBytes(const pybind11::bytes& bytes) {
  return Scine::molassembler::JSONSerialization(
    binaryFromPybind11Bytes(bytes),
    Scine::molassembler::JSONSerialization::BinaryFormat::CBOR
  );
}

Scine::molassembler::Molecule fromBSONBytes(const pybind11::bytes& bytes) {
  return Scine::molassembler::JSONSerialization(
    binaryFromPybind11Bytes(bytes),
    Scine::molassembler::JSONSerialization::BinaryFormat::BSON
  );
}

void init_serialization(pybind11::module& m) {
  using namespace Scine::molassembler;

  m.def(
    "to_json",
    &toJSON,
    pybind11::arg("molecule"),
    "Serialize a molecule into a JSON string"
  );
  m.def(
    "to_cbor",
    &toCBORBytes,
    pybind11::arg("molecule"),
    "Serialize a molecule into CBOR JSON"
  );
  m.def(
    "to_bson",
    &toBSONBytes,
    pybind11::arg("molecule"),
    "Serialize a molecule into BSON JSON"
  );
  m.def(
    "from_json",
    &fromJSON,
    pybind11::arg("json_str"),
    "Deserialize a JSON string into a molecule instance"
  );
  m.def(
    "from_cbor",
    &fromCBORBytes,
    pybind11::arg("cbor_bytes"),
    "Deserialize CBOR bytes into a molecule instance"
  );
  m.def(
    "from_bson",
    &fromBSONBytes,
    pybind11::arg("bson_bytes"),
    "Deserialize BSON bytes into a molecule instance"
  );
}
