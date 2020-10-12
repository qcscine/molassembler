__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

import os
from conans import ConanFile, CMake, tools


def conan_paths_str(build_folder):
    return os.path.join(build_folder, "conan_paths.cmake")


class TestPackageConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake_paths"
    exports_sources = "CMakeLists.txt", "test.cpp"

    _cmake = None

    def build_requirements(self):
        if not tools.which("cmake") or CMake.get_version() < "3.13.4":
            self.build_requires("cmake/[>3.13.4]@scine/stable")

    def _configure(self):
        if self._cmake:
            return self._cmake

        self._cmake = CMake(self)
        self._cmake.definitions["CMAKE_PROJECT_MolassemblerTestPackage_INCLUDE"] = conan_paths_str(
            self.build_folder)
        self._cmake.configure()
        return self._cmake

    def build(self):
        cmake = self._configure()
        cmake.build()

    def test(self):
        cmake = self._configure()
        cmake.test(output_on_failure=True)
