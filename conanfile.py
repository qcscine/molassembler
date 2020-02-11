from conans import ConanFile, CMake, tools
from pathlib import Path
from textwrap import wrap
import sys


def conan_paths(build_folder):
    """Gets the full file path of the generated conan paths"""
    return str(Path(build_folder) / "conan_paths.cmake")


def python_module_dir(package_folder):
    """Generates the package's installed python path"""
    python_dir = "python" + str(sys.version_info.major) + \
        "." + str(sys.version_info.minor)
    return str(Path(package_folder) / "lib" / python_dir / "site-packages")


class MolassemblerConan(ConanFile):
    name = "molassembler"
    version = "1.0.0"
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
    options = {"shared": [True, False], "python": [True, False]}
    default_options = {"shared": True, "python": False}
    generators = "cmake_paths"
    exports_sources = "src/*", "CMakeLists.txt"
    requires = [("boost/[~=1.71.0]@conan/stable"),
                ("scine_utilities/[~=2.1.0]"),
                ("eigen/[~=3.3.7]@conan/stable")]

    def _configure_cmake(self):
        cmake = CMake(self)
        additional_definitions = {
            "SCINE_MARCH": "",
            "SCINE_BUILD_DOCS": False,
            "SCINE_BUILD_TESTS": False,
            "SCINE_BUILD_PYTHON_BINDINGS": self.options.python,
            "PYTHON_EXECUTABLE": sys.executable,
            "CMAKE_PROJECT_molassembler_INCLUDE": conan_paths(self.build_folder)
        }
        cmake.definitions.update(additional_definitions)
        # Mess with cmake definitions, etc.
        cmake.configure()
        return cmake

    def configure(self):
        tools.check_min_cppstd(self, "14")

    def build_requirements(self):
        if self.settings.os == "Windows":
            self.build_requires("mingw_installer/[~=1.0]@conan/stable")

        # TODO: As soon as 2.4.2 is available on conan, prefer this
        # if self.options.python:
        #     self.build_requires("pybind11/2.4.2@conan/stable")

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()
        cmake.patch_config_paths()

    def package_info(self):
        self.cpp_info.libs = ["molassembler"]
        if self.options.python:
            self.env_info.PYTHONPATH.append(
                python_module_dir(self.package_folder)
            )
