=======
Scripts
=======

This script is provided to ease installations in offline environments such as
scientific computing clusters. Molassembler's CMake code pulls some files from
the internet during configuration. This can be avoided by specifying the CMake
option ``-DMOLASSEMBLER_OFFLINE_RESOURCE_DIR=/absolute/path/to/resources``.

Usage::

    $ bash download-resources.sh
    $ scp -r resources user@target-machine:/where/to/store/resources
    $ ssh user@target-machine
    $ cd /path/to/molassembler/build
    $ cmake -DMOLASSEMBLER_OFFLINE_RESOURCE_DIR=/where/to/store/resources ..
