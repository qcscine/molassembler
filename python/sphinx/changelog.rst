Changelog
=========

This project adheres to semantic versioning.

This changelog only captures changes to the python bindings. For detailed
changes to the underlying C++ library, see the repository's ``CHANGELOG.rst``.

1.2.1
-----

No Python-specific changes

1.2.0
-----

No Python-specific changes

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
