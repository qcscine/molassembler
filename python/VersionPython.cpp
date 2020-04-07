/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "molassembler/Version.h"

void init_version(pybind11::module& m) {
  using namespace Scine::molassembler;

  m.attr("__version__") = Scine::molassembler::version::fullVersion();

  auto versionSubmodule = m.def_submodule("version");

  versionSubmodule.attr("MAJOR") = pybind11::int_(version::major);
  versionSubmodule.attr("MINOR") = pybind11::int_(version::minor);
  versionSubmodule.attr("PATCH") = pybind11::int_(version::patch);

  versionSubmodule.def(
    "major_minor",
    &version::majorMinor,
    "Returns a major.minor formatted string of the molassembler version"
  );

  versionSubmodule.def(
    "full_version",
    &version::fullVersion,
    "Returns a major.minor.patch formatted string of the molassembler version"
  );

  versionSubmodule.def(
    "compiled",
    []() -> std::string {
      return std::string(__DATE__) + " " + std::string(__TIME__);
    },
    "Returns a string of date and time the python bindings were compiled"
  );
}
