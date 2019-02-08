import setuptools

# Read README.md for the long description
with open("README.md", "r") as fh:
  long_description = fh.read()

# Define the setup
setuptools.setup(
  name="molassembler",
  version="0.1.0",
  author="Jan-Grimo Sobez",
  author_email="jan-grimo.sobez@phys.chem.ethz.ch",
  description="Inorganic molecular graph interpretation, modification and conformer generation",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://gitlab.chab.ethz.ch/scine/molassembler",
  packages=[""],
  package_dir={"": "."},
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
