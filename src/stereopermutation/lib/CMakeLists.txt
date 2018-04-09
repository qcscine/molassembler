set(stereopermutation_VERSION 0.1.0)

# dependencies
find_package(Eigen3 3.1.1 REQUIRED)
find_package(chemical_symmetries 0.1.0 REQUIRED)

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Stereopermutation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/GenerateUniques.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/LogicalOperatorTests.h
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Util.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/Stereopermutation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation/GenerateUniques.cpp
)

add_library(stereopermutation STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(stereopermutation PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(stereopermutation SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)
target_include_directories(stereopermutation PUBLIC ${EIGEN3_INCLUDE_DIRS})
target_link_libraries(stereopermutation chemicalsymmetries)

if(UNIX)
  add_definitions(-fPIC)
endif()

target_compile_options(stereopermutation PUBLIC -Wall -Wpedantic -Wextra)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  message(STATUS "Selected release build")
  # no debug information, disable assert
  target_compile_definitions(stereopermutation PRIVATE NDEBUG)
  target_compile_options(stereopermutation PRIVATE -O2)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(stereopermutation PRIVATE -g -O2)
endif()

# Install all relevant header files in <installdir>/include
include(CMakePackageConfigHelpers)

# Set destination for the generated CMake config files
set(ConfigPackageLocation lib/cmake/stereopermutation)
set(INCLUDE_INSTALL_DIR include)

# Generate a basic configuration file in the binary director and subdir
# stereopermutation
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/stereopermutation/stereopermutationConfigVersion.cmake"
  VERSION ${stereopermutation_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/stereopermutation/stereopermutationConfig.cmake"
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS INCLUDE_INSTALL_DIR
)

install(
  TARGETS stereopermutation
  EXPORT stereopermutationTargets
  DESTINATION lib
)

install(
  EXPORT stereopermutationTargets
  DESTINATION "${ConfigPackageLocation}"
)

install(DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/stereopermutation"
  DESTINATION
  ${INCLUDE_INSTALL_DIR}
  FILES_MATCHING
  PATTERN "*.h"
)

# Install ...Config.cmake and ...ConfigVersion into specified directory
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/stereopermutation/stereopermutationConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/stereopermutation/stereopermutationConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
)

# Source groups for Visual Studio
source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
