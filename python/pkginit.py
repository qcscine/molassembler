# This is more or less copy-paste from importlib's docs, makes the pybind11
# module .so wild-importable (e.g. from `scine_molassembler.io import *`)
import importlib.util
import os
import sys
from pathlib import Path
spec = importlib.util.spec_from_file_location(
    __name__,
    os.path.join(Path(__file__).parent.absolute(), __name__ + ".so")
)
module = importlib.util.module_from_spec(spec)
sys.modules[__name__] = module
spec.loader.exec_module(module)
