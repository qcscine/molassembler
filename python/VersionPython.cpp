/*!@file
 * @copyright This code is licensed under the 3-clause BSD license.
 *   Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *   See LICENSE.txt for details.
 */
#include "TypeCasters.h"
#include "Molassembler/Version.h"

void init_version(pybind11::module& m) {
  using namespace Scine::Molassembler;

  m.attr("__version__") = Scine::Molassembler::Version::full();

  auto versionSubmodule = m.def_submodule("version");

  versionSubmodule.attr("MAJOR") = pybind11::int_(Version::major);
  versionSubmodule.attr("MINOR") = pybind11::int_(Version::minor);
  versionSubmodule.attr("PATCH") = pybind11::int_(Version::patch);

  versionSubmodule.def(
    "major_minor",
    &Version::majorMinor,
    "Returns a major.minor formatted string of the molassembler version"
  );

  versionSubmodule.def(
    "full_version",
    &Version::full,
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
