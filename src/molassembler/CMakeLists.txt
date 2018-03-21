# oldest tested version of cmake used, feel free to lower this if older versions
# work too
cmake_minimum_required(VERSION 3.5.1)

project(molassembler)

# project declaration
enable_testing()

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(analysis)
