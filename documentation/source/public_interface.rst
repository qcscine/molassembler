Public interface
================

v1.0.0
------
All headers in the base library directory of molassembler and files included
from there form the public interface. None of the files in deeper directories
are part of the public interface.

The auxiliary library headers that are by inclusion also part of the public
interface are currently:

- chemical symmetries/

  - Names

- temple/

  - Preprocessor
  - TinySet
  - Traits
  - constexpr/

    - Bitmask
    - Math

It is inadvisable to include these files directly and make use of their
functionality since the types therein may be moved around in future versions.
Any molassembler headers that use types from these files include them, and there
are no forward declarations.
