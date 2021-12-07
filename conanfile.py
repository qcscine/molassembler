__copyright__ = """This code is licensed under the 3-clause BSD license.
Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
See LICENSE.txt for details.
"""

from dev.conan.base import ScineConan


class MolassemblerConan(ScineConan):
    name = "scine_molassembler"
    version = "1.0.0"
    license = "BSD-3-Clause"
    author = "Jan-Grimo Sobez jan-grimo.sobez@phys.chem.ethz.ch"
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
        "microarch": "none"
    }
    exports = "dev/conan/*.py"
    exports_sources = [
        "dev/conan/hook.cmake",
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
        "scine_utilities/[~=3.0.0]@scine/develop",
        "RingDecomposerLib/1.1.3@scine/stable",
        "nauty/2.7r1@scine/stable"
    ]

    def _configure_cmake(self):
        return super()._configure_cmake_base("molassembler", None)
