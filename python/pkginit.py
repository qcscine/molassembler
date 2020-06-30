import importlib.util
import os
import sys
from pathlib import Path
from distutils import sysconfig

expected_suffix = sysconfig.get_config_var("EXT_SUFFIX")
expected_name = __name__ + expected_suffix
this_file_dir = Path(__file__).parent.absolute()
python_module_path = os.path.join(this_file_dir, expected_name)
if not os.path.exists(python_module_path):
    raise ImportError("Could not find {}".format(expected_name))

# This is more or less copy-paste from importlib's docs, makes the pybind11
# module wild-importable (e.g. from `scine_molassembler.IO import *`)
spec = importlib.util.spec_from_file_location(__name__, python_module_path)
module = importlib.util.module_from_spec(spec)
sys.modules[__name__] = module
spec.loader.exec_module(module)
