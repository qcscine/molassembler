/*!@file
 * @copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt
 */
#include "pybind11/pybind11.h"
#include "molassembler/Version.h"

void init_version(pybind11::module& m) {
  using namespace Scine::molassembler;

  auto versionSubmodule = m.def_submodule("version");

  versionSubmodule.attr("MAJOR") = pybind11::int_(version::major);
  versionSubmodule.attr("MINOR") = pybind11::int_(version::minor);
  versionSubmodule.attr("PATCH") = pybind11::int_(version::patch);

  versionSubmodule.def(
    "major_minor",
    &version::majorMinor
  );

  versionSubmodule.def(
    "full_version",
    &version::fullVersion
  );

  versionSubmodule.def(
    "compiled",
    []() -> std::string {
      return std::string(__DATE__) + " " + std::string(__TIME__);
    }
  );
}
