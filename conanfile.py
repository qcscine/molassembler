import os
import sys
import re
import subprocess as sp
from conans import ConanFile, CMake, tools
from conans.errors import ConanInvalidConfiguration


def find_arch_in_output(cmdlist, regex):
    result = sp.run(cmdlist, stdout=sp.PIPE,
                    stderr=sp.STDOUT, universal_newlines=True)
    result.check_returncode()
    matcher = re.compile(regex)

    for match in matcher.finditer(result.stdout):
        return match.group("arch")

    for match in matcher.finditer(result.stderr):
        return match.group("arch")

    return None


def determine_arch(conanfile):
    arch = None

    if conanfile.settings.compiler == "gcc":
        arch = find_arch_in_output(
            ["gcc", "-march=native", "-Q", "--help=target"],
            r"-march=\s+(?P<arch>[A-z0-9]+)"
        )

    if conanfile.settings.compiler in ["clang", "apple-clang"]:
        arch = find_arch_in_output(
            ["clang", "-march=native", "-xc", "-", "-###"],
            "\"-target-cpu\"\\s+\"(?P<arch>[A-z0-9]+)\""
        )

    return arch or ""


def cmake_at_least(conanfile, version):
    if not tools.which("cmake"):
        return False

    return CMake.get_version() >= version


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
    description = "Molecular graph interpretation, modification and conformer generation"
    topics = ("chemistry", "cheminformatics", "molecule")
    settings = "os", "compiler", "build_type", "arch"
    options = {
        "shared": [True, False],
        "python": [True, False],
        "tests": [True, False],
        "docs": [True, False],
        "coverage": [True, False],
        "microarch": ["detect", "none"]
    }
    default_options = {
        "shared": True,
        "python": False,
        "tests": False,
        "docs": False,
        "coverage": False,
        "microarch": "none"
    }
    generators = "cmake"
    exports_sources = [
        "cmake/*",
        "src/*",
        "doc/*",
        "test/*",
        "python/*",
        "CMakeLists.txt",
        ".conan_include.cmake",
    ]
    requires = "scine_utilities/[~=2.1.0]@ci/develop"
    revision_mode = "scm"

    _cmake = None
    _microarch = ""

    def _ensure_microarch_set(self):
        if not self._microarch and self.options.microarch == "detect":
            self._microarch = determine_arch(self)
            if self._microarch != "":
                self.output.info(
                    "Detected microarch: {}".format(self._microarch))

    def _configure_cmake(self):
        if self._cmake:
            return self._cmake

        self._ensure_microarch_set()

        self._cmake = CMake(self)
        additional_definitions = {
            "SCINE_MARCH": self._microarch,
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

        if self.options.microarch == "none":
            self.options["scine_utilities"].microarch = "none"

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

        self._ensure_microarch_set()
        self.info.options["microarch"] = self._microarch

    def build_requirements(self):
        if self.options.python:
            self.build_requires("pybind11/2.4.2@scine/dependencies")

        if not cmake_at_least(self, "3.13.4"):
            self.build_requires("cmake_installer/[~=3.13.4]@conan/stable")

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

        if self.options.tests:
            # Run tests sequentially
            cmake.parallel = False
            cmake.test(output_on_failure=True)
            cmake.parallel = True

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
