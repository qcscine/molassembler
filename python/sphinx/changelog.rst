Changelog
=========

This project adheres to semantic versioning.

This changelog only captures changes to the python bindings. For detailed
changes to the underlying C++ library, see the repository's ``CHANGELOG.rst``.

1.1.0
-----
Added
~~~~~
- Added modifying functions to ``Graph``
- Direct copying support for ``Molecule`` instead of via pickling
- Added build-time generation of typing stubs with pybind11-stubgen

Changed
~~~~~~~
- Altered name of ``ChiralStatePreservation`` enum member from ``None`` to
  ``DoNotPreserve`` (the former is a reserved keyword)
- Better automatic type signature annotations in docstrings

1.0.0 - 2020-01-23
------------------
Added
~~~~~
- Doctested examples
- Pickle support for ``Molecule`` using serialization
- Experimental library-internal SMILES parser in ``io.experimental``
- Molecule instances integrate nicely with notebooks using ``_repr_svg_``
- ``__repr__`` members for ``Molecule``, ``Graph`` and ``StereopermutatorList``
- Improved ``__init__.py`` to allow wild imports
- ``setup.py`` correctly installs as platlib

Changed
~~~~~~~
- The ``interpret`` family of functions was renamed ``molecules`` and moved to a
  new submodule ``interpret`` in lieu of the addition of very similar functions
  that yield only ``Graph`` instances without constructing ``Molecule``
  instances.
- The following member functions of ``AtomStereopermutator`` and
  ``BondStereopermutator`` are now properties:

  - ``assigned``
  - ``index_of_permutation``
  - ``num_assignments``
  - ``num_stereopermutations``

- The member functions ``N`` and ``B`` of ``Graph`` are now properties
