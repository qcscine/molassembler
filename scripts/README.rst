=======
Scripts
=======

These scripts are provided as convenience CMake installation helpers for the
RingDecomposerLib and nauty libraries. Both of these libraries' distributions
do not currently have CMake wrapper code suitable for CMake's ``find_package``
functionality, so both libraries are intrusively modeled or modified to provide
configuration files.

Molassembler pulls the source files for both libraries during
installation if they are not available via ``find_package``. For installations
without internet access (e.g. scientific computing clusters for security
reasons), it is therefore advisable to use these scripts to install these
dynamic dependencies and point to their installation directories with
``-DCMAKE_PREFIX_PATH`` when installing Molassembler itself.


How to use
==========

1. Run each script once without arguments on a machine with an internet
   connection. They will pull their dependencies and terminate.
2. Copy the archive, the installation script, and any name-related patch or
   cmake files to the target machine.
3. Run the install script, supplying the installation prefix as an extra
   unflagged argument

For instance, for nauty::

    $ bash install-nauty.sh
    $ scp install-nauty.sh nauty.cmake nauty27r1.tar.gz user@target-machine:/path
    $ ssh user@target-machine
    $ cd /path
    $ bash install-nauty.sh /path/to/install/into
