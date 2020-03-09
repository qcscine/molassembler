import os
import sys
from conans import ConanFile, CMake, tools
from conans.errors import ConanInvalidConfiguration


def python_module_dir(package_folder):
    """Generates the package's installed python path"""
    python_dir = "python" + str(sys.version_info.major) + \
        "." + str(sys.version_info.minor)
    return os.path.join(package_folder), "lib", python_dir, "site-packages"


class MolassemblerConan(ConanFile):
    name = "scine_molassembler"
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
    options = {key: [True, False]
               for key in ["shared", "python", "tests", "docs", "coverage"]}
    default_options = {
        "shared": True,
        "python": False,
        "tests": False,
        "docs": False,
        "coverage": False
    }
    generators = "cmake"
    exports_sources = "src/*", "CMakeLists.txt", ".conan_include.cmake", "doc/*"
    build_requires = [("cmake_installer/[~=3.13.4]@conan/stable")]
    requires = [("scine_utilities/[~=2.1.0]@ci/develop")]
    revision_mode = "scm"

    _cmake = None

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake

        self._cmake = CMake(self)
        additional_definitions = {
            "SCINE_MARCH": "",
            "SCINE_BUILD_DOCS": self.options.docs,
            "SCINE_BUILD_TESTS": self.options.tests,
            "SCINE_BUILD_PYTHON_BINDINGS": self.options.python,
            "PYTHON_EXECUTABLE": sys.executable,
            "CMAKE_PROJECT_molassembler_INCLUDE": ".conan_include.cmake",
            "COVERAGE": self.options.coverage
        }
        self._cmake.definitions.update(additional_definitions)
        # Mess with cmake definitions, etc.
        self._cmake.configure()
        return self._cmake

    def configure(self):
        tools.check_min_cppstd(self, "14")

        if self.options.python:
            self.options["scine_utilities"].python = True

        if self.options.coverage:
            if not self.options.tests:
                raise ConanInvalidConfiguration(
                    "Coverage testing requires testing to be enabled"
                )

            if self.settings.build_type != "Debug":
                raise ConanInvalidConfiguration(
                    "Coverage testing should be done on a debug build")

    def package_id(self):
        # Tests and docs do not contribute to package id
        del self.info.options.tests
        del self.info.options.docs
        del self.info.options.coverage

    def build_requirements(self):
        if self.options.python:
            self.build_requires("pybind11/2.4.2@scine/dependencies")

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

        if self.options.tests:
            cmake.test(output_on_failure=True)

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()
        cmake.patch_config_paths()

    def package_info(self):
        self.cpp_info.libs = tools.collect_libs(self)
        if self.options.python:
            self.env_info.PYTHONPATH.append(
                python_module_dir(self.package_folder)
            )
