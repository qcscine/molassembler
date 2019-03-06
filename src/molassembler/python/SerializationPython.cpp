/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "molassembler/Molecule.h"
#include "molassembler/Serialization.h"

void init_serialization(pybind11::module& m) {
  using namespace Scine::molassembler;

  m.def(
    "to_json",
    &toJSON,
    pybind11::arg("molecule"),
    "Serialize a molecule into a JSON string"
  );
  m.def(
    "to_base64",
    &toBase64EncodedCBOR,
    pybind11::arg("molecule"),
    "Serialize a molecule into a base 64 string"
  );
  m.def(
    "from_json",
    &fromJSON,
    pybind11::arg("json_str"),
    "Deserialize a JSON string into a molecule instance"
  );
  m.def(
    "from_base64",
    &fromJSON,
    pybind11::arg("base64_str"),
    "Deserialize a base 64 string into a molecule instance"
  );
}
