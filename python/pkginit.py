def bootstrap():
    # This is more or less copy-paste from importlib's docs
    import importlib.util
    import os
    import sys
    from pathlib import Path

    path = os.path.join(
        Path(__file__).parent.absolute(),
        "scine_molassembler.so"
    )
    module_name = "scine_molassembler"

    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)

    # Erase the function itself
    global bootstrap
    del bootstrap


bootstrap()
