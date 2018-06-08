# TODO
- Clean up CMake: Reduce code duplication -> see all the installs
- Separate constexpr temple tests just like the non-constexpr ones
- Before publishing this in any way, must consider visibility of core functions
  and/or making only a considerably more limited header file set public
  (i.e. only Molecule, generateConformation, StereocenterList, Stereocenter,
  CNStereocenter and EZStereocenter public function headers)
  see particularly: https://gcc.gnu.org/wiki/Visibility
- Have a look at link time optimization and profile guided optimization so that
  the parts constituting libmolassembler.a fuse nicely
