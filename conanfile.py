from conans import ConanFile, CMake
from pathlib import Path
from textwrap import wrap


def conan_paths_path_str(build_folder):
    return str(Path(build_folder) / "conan_paths.cmake")


class MolassemblerConan(ConanFile):
    name = "Molassembler"
    version = "0.1.0"
    license = "BSD-3-Clause"
    author = "Jan-Grimo Sobez jan-grimo.sobez@phys.chem.ethz.ch"
    url = "https://gitlab.chab.ethz.ch/scine/molassembler"
    description = """
Molassembler is a C++ library that aims to facilitate crossings between
Cartesian and graph representations of molecules. It provides the necessary
functionality to represent a molecule as a graph, modify it in graph space, and
generate coordinates from graphs. It can capture the absolute configuration
of multidentate and haptic inorganic molecules from positional data and
generate non-superposable stereopermutations as output."""
    topics = ("chemistry", "cheminformatics", "molecule")
    settings = "os", "compiler", "build_type", "arch"
    options = {"subsume_dependencies": [True, False]}
    default_options = {"subsume_dependencies": False}
    generators = "cmake_paths"
    exports_sources = "src/*", "CMakeLists.txt"
    build_requires = [("eigen/[~=3.3.7]@conan/stable")]
    requires = [("boost/[~=1.71.0]@conan/stable"),
                ("lapack/[~=3.7.1]@conan/stable")]

    def _configure_cmake(self):
        cmake = CMake(self)
        additional_definitions = {
            "SCINE_BUILD_DOCS": False,
            "SCINE_BUILD_TESTS": False,
            "SCINE_BUILD_PYTHON_BINDINGS": False,
            "CMAKE_PROJECT_molassembler_INCLUDE": conan_paths_path_str(self.build_folder)
        }
        cmake.definitions.update(additional_definitions)
        # Mess with cmake definitions, etc.
        cmake.configure()
        return cmake

    def build_requirements(self):
        if self.settings.os == "Windows":
            self.build_requires("mingw_installer/[~=1.0]@conan/stable")

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()
        cmake.patch_config_paths()

    def package_info(self):
        self.cpp_info.libs = ["molassembler"]


if __name__ == "__main__":
    print("Molassembler conan recipe options and default values:\n")

    explanations = {
        "subsume_dependencies": "Use static dependencies in order to get libraries and binaries with minimal runtime dependencies."
    }

    recipe = MolassemblerConan
    for key in recipe.options.keys():
        print("- " + str(key) + ": " + str(recipe.default_options[key]))
        if key in explanations.keys():
            for line in wrap(explanations[key], initial_indent="  ", subsequent_indent="  "):
                print(line)

            print()