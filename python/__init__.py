def __bootstrap__():
    global __bootstrap__, __loader__, __file__
    import pkg_resources as pkg
    import imp
    __file__ = pkg.resource_filename(__name__, 'scine_molassembler.so')
    __loader__ = None
    del __bootstrap__, __loader__
    imp.load_dynamic(__name__, __file__)


__bootstrap__()
