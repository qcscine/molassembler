# Contribution guidelines

## Code style

Molassembler follows SCINE code style conventions.

### signed vs unsigned
Molassembler prefers unsigned integers for container element indices (follow
STL style).

- Generally prefer range-for over indexed access
- Iterator differences return signed types

Ideally, container element indices should be `std::size_t` and its iterator
differences should be `std::ptrdiff_t` if they are expected to be able to
contain large amounts of values.

Any integers where subtracting is expected should be signed. If unsigned values
are used in subtractions, `static_cast` them beforehand!

### C++ style
Try to follow the C++ Core Guidelines.

### Formatting
Unenforced, follows my own idiosyncratic formatting style that I think
clang-format doesn't support. I try to keep the code 80-char wide.
