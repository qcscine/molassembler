import setuptools

# Read README.md for the long description
markdown_description = """
# Molassembler library
## Overview

Molassembler is a C++ library that aims to facilitate crossings between
Cartesian and graph representations of molecules. It provides the necessary
functionality to represent a molecule as a graph, modify it in graph space, and
generate coordinates from graphs. It can capture the absolute configuration
of multidentate and haptic inorganic molecules from positional data and
generate non-superposable stereopermutations as output.
"""

# Define the setup
setuptools.setup(
    name="molassembler",
    version="0.1.0",
    author="Jan-Grimo Sobez",
    author_email="jan-grimo.sobez@phys.chem.ethz.ch",
    description="Inorganic molecular graph interpretation, modification and conformer generation",
    long_description=markdown_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.chab.ethz.ch/scine/molassembler",
    packages=[""],
    package_data={"": ["molassembler.so"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: Other/Proprietary License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    zip_safe=False,
    test_suite='pytest',
    tests_require=['pytest']
)
