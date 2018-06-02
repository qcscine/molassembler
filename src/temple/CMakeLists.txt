# How to make this into a library?

# oldest tested version of cmake used (too lazy to figure out newest compatible
# cmake version)
cmake_minimum_required(VERSION 3.5.1)

# declaration
project(temple VERSION 0.1.0)

# Add testing
enable_testing()

# Subdirectories
add_subdirectory(lib)
add_subdirectory(test)
