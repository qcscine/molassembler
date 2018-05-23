#ifndef INCLUDE_MOLASSEMBLER_SERIALIZE_H
#define INCLUDE_MOLASSEMBLER_SERIALIZE_H

#include "Molecule.h"

namespace molassembler {

std::string toJSON(const Molecule& molecule);
std::vector<std::uint8_t> toCBOR(const Molecule& molecule);
std::string toBase64EncodedCBOR(const Molecule& molecule);

//! Deserialize a JSON string representation of a Molecule
Molecule fromJSON(const std::string& serializedMolecule);
Molecule fromCBOR(const std::vector<std::uint8_t>& cbor);
Molecule fromBase64EncodedCBOR(const std::string& base64EncodedCBOR);

} // namespace molassembler

#endif
