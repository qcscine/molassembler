# Contribution guidelines

## Structure
Molassembler widely uses the pImpl pattern.

Headers in the top-level library directory are meant for library consumers. Any
subfolder headers are implementation details and do not form part of the public
interface. No top-level headers may therefore include sublevel headers.


## Code style

Molassembler follows SCINE code style convention.

### signed vs unsigned
Molassembler prefers using unsigned to represent container element indices
(follow STL style).

- Generally prefer range-for over indexed access
- Iterator differences return signed types
- Do not waste bit representation space

Ideally, container element indices should be `std::size_t` and its iterator
differences should be `std::ptrdiff_t` if they are expected to be able to
contain large amounts of values.

Any integers where subtracting is expected should be signed. If unsigned values
are used in subtractions, `static_cast` them!

### C++ style
Try to follow the C++ Core Guidelines.

### Formatting
Unenforced since I currently have my own idiosyncratic formatting style that I
think clang-format doesn't support. I try to keep the code 80-char wide at most.
