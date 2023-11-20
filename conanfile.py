__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
See LICENSE.txt for details.
"""

from dev.conan.base import ScineConan


class MolassemblerConan(ScineConan):
    name = "scine_molassembler"
    version = "2.0.1"
    license = "BSD-3-Clause"
    author = "ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group"
    url = "https://scine.ethz.ch/download/molassembler"
    description = "Molecular graph interpretation, modification and conformer generation"
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
        "microarch": "none",
        "RingDecomposerLib:shared": False,
        "nauty:shared": False
    }
    exports = "dev/conan/*.py"
    exports_sources = [
        "dev/conan/hook.cmake",
        "dev/conan/glue/*",
        "dev/cmake/*",
        "cmake/*",
        "extern/*",
        "src/*",
        "doc/*",
        "test/*",
        "python/*",
        "CMakeLists.txt",
    ]
    requires = [
        "scine_utilities/9.0.0",
        "RingDecomposerLib/1.1.3",
        "nauty/2.7r1"
    ]
    cmake_name = "Molassembler"

    def package_info(self):
        super().package_info()

        self.cpp_info.components["Molassembler"].cxxflags = ["-fopenmp"]
        self.cpp_info.components["Molassembler"].sharedlinkflags = ["-fopenmp"]
        self.cpp_info.components["Molassembler"].exelinkflags = ["-fopenmp"]
