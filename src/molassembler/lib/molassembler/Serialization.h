/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 * @brief Provides serialization / deserialization for Molecule instances
 */

#ifndef INCLUDE_MOLASSEMBLER_SERIALIZATION_H
#define INCLUDE_MOLASSEMBLER_SERIALIZATION_H

#include <vector>
#include <string>

namespace Scine {

namespace molassembler {

// Forward-declarations
class Molecule;

/*!
 * @brief Serializes a molecule into a compact JSON representation
 *
 * The JSON representation, although principally human readable, is very
 * compact. Keys are heavily shortened. Each Molecule JSON Object has the
 * following structure:
 *
 * - a: List of AtomStereopermutator Objects
 *   - a: Assignment index (key omitted if unassigned)
 *   - c: Central index
 *   - r: Ranking
 *     - s: Sorted subtituents
 *     - l: Ligands
 *     - lr: Ranked ligands
 *     - lnk: Links (key omitted if empty)
 *   - s: Symmetry index
 * - b: List of BondStereopermutator Objects
 *   - a: Assignment index (key omitted if unassigned)
 *   - e: Edge on which it is placed
 * - g: Graph Object
 *   - Z: List of atomic numbers
 *   - E: List of edges, each a List
 *     - 0: Source vertex
 *     - 1: Target vertex
 *     - 2: Bond type index
 * - v: Library version List
 *   - 0: Major
 *   - 1: Minor
 *   - 2: Patch
 *
 * @returns An unprettified JSON string serialization
 */
std::string toJSON(const Molecule& molecule);
/*!
 * @brief Serializes a molecule into compact binary object notation
 * @returns Variable-length binary
 */
std::vector<std::uint8_t> toCBOR(const Molecule& molecule);

/*!
 * @brief Serializes a molecule into base 64 encoded compact binary object
 *   representation
 */
std::string toBase64EncodedCBOR(const Molecule& molecule);

//! Deserialize a JSON string representation of a Molecule
Molecule fromJSON(const std::string& serializedMolecule);
//! Deserialize a compact binary object notation serialization of a molecule
Molecule fromCBOR(const std::vector<std::uint8_t>& cbor);
//! Deserialize a base 64 encoded compact binary object representation of a molecule
Molecule fromBase64EncodedCBOR(const std::string& base64EncodedCBOR);

} // namespace molassembler

} // namespace Scine

#endif
