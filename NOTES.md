# TODO

- Clean up CMake: Reduce code duplication -> see all the installs
- Separate constexpr temple tests just like the non-constexpr ones
- Figure out how to integrate all documentations with each other
- Isn't it testable to see if the tetrahedron definitions really encompass all
  chiral elements in the symmetry? I think this is tested in the library
  itself, too
- Before publishing this in any way, must consider visibility of core functions
  and/or making only a considerably more limited header file set public
  (i.e. only Molecule, generateConformation, StereocenterList, Stereocenter,
  CNStereocenter and EZStereocenter public function headers)
  see particularly: https://gcc.gnu.org/wiki/Visibility
- Add a GOR1 test!
