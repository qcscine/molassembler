Contributing to SCINE Molassembler
==================================

Contribution Process
--------------------

The development for this code is done in a private repository maintained by the
Reiher Research Group. GitHub is only used for the official releases.

If you would like to contribute a larger change, please write to scine@phys.chem.ethz.ch.
For smaller changes, you can create a pull request on GitHub. If we agree with
the changes, a member of the Reiher Research Group will include them in our
development code. Of course, we will give proper acknowledgment for any external
contribution (see below for a list of all contributors). As soon as these changes
are available in an official release, we will close the corresponding pull requests
and/or issues on GitHub.

Please note that contributing a small change does in no way mean that you will
be added to the author list of a future paper and/or Zenodo entry!

Main Contributors
-----------------

Almost all contributions to SCINE in general and this repository in specific come
from members of the Reiher research group.

Further Contributors
--------------------

So far, no one else has contributed to this repository

Contribution Guidelines
-----------------------

Code style
..........

Molassembler follows SCINE code style conventions.

signed vs unsigned
..................

Molassembler prefers unsigned integers for container element indices (follow
STL style).

- Generally prefer range-for over indexed access
- Iterator differences return signed types

Ideally, container element indices should be ``std::size_t`` and its iterator
differences should be ``std::ptrdiff_t`` if they are expected to be able to
contain large amounts of values.

Any integers where subtracting is expected should be signed. If unsigned values
are used in subtractions, ``static_cast`` them beforehand!

C++ style
.........

Try to follow the C++ Core Guidelines.

Formatting
..........

Unenforced, follows my own idiosyncratic formatting style that I think
clang-format doesn't support. I try to keep the code 80-char wide.
