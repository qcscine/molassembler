import os
from conans import ConanFile, CMake, tools


def conan_paths_str(build_folder):
    return os.path.join(build_folder, "conan_paths.cmake")


class TestPackageConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake_paths"
    exports_sources = "CMakeLists.txt", "test.cpp"

    def build_requirements(self):
        if not tools.which("cmake") or CMake.get_version() < "3.13.4":
            self.build_requires("cmake_installer/[~=3.13.4]@conan/stable")

    def _configure(self):
        cmake = CMake(self)
        cmake.definitions["CMAKE_PROJECT_MolassemblerTestPackage_INCLUDE"] = conan_paths_str(
            self.build_folder)
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure()
        cmake.build()

    def test(self):
        cmake = self._configure()
        cmake.test()
