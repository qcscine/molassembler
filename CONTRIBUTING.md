# Contribution guidelines

## Structure
Any headers (and their corresponding implementations) that are meant for
use by library consumers MUST be placed into the top-level library directory.
Any headers that are NOT meant for consumer access but are auxiliary to
top-level implementations MUST be in a subfolder.

No top-level header may include a subfolder header.


## Code style

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
Try to follow the C++ Core Guidelines. Clang-tidy will enforce some aspects of
it.

### Formatting
Enforced by clang-format.
