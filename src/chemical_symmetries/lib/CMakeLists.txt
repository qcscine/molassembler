set(chemical_symmetries_VERSION 0.1.0)

# dependencies
find_package(Eigen3 3.1.1 REQUIRED)
find_package(temple 0.1.0 REQUIRED)

message(STATUS "temple include dirs: ${TEMPLE_INCLUDE_DIR}")

set(HEADERS_All
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/AngleLookup.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/CompileTimeOptions.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/ConstexprProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/DynamicProperties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Primitives.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Properties.h
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Symmetries.h
)

set(SOURCES_All
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/DynamicProperties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Primitives.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries/Symmetries.cpp
)

add_library(chemicalsymmetries STATIC ${SOURCES_All} ${HEADERS_all})
target_include_directories(
  chemicalsymmetries
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(
  chemicalsymmetries
    SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
target_include_directories(
  chemicalsymmetries
    PUBLIC
      ${EIGEN3_INCLUDE_DIRS}
      ${TEMPLE_INCLUDE_DIR}
)
target_link_libraries(chemicalsymmetries temple)

if(UNIX)
  add_definitions(-fPIC)
endif()

target_compile_options(chemicalsymmetries PUBLIC -Wall -Wpedantic -Wextra)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(
    chemicalsymmetries
    PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0
  )
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  message(STATUS "Selected release build")
  # no debug information, disable assert
  target_compile_definitions(chemicalsymmetries PRIVATE NDEBUG)
  target_compile_options(chemicalsymmetries PRIVATE -O2)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(chemicalsymmetries PRIVATE -g -O2)
endif()

# Install all relevant header files in <installdir>/include
include(CMakePackageConfigHelpers)

# Set destination for the generated CMake config files
set(ConfigPackageLocation lib/cmake/chemical_symmetries)
set(INCLUDE_INSTALL_DIR include)

# Generate a basic configuration file in the binary director and subdir
# chemical_symmetries
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/chemical_symmetries/chemical_symmetriesConfigVersion.cmake"
  VERSION ${chemical_symmetries_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/chemical_symmetries/chemical_symmetriesConfig.cmake"
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS INCLUDE_INSTALL_DIR
)

install(
  TARGETS chemicalsymmetries
  EXPORT chemical_symmetriesTargets
  DESTINATION lib
)

install(
  EXPORT chemical_symmetriesTargets
  DESTINATION "${ConfigPackageLocation}"
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/chemical_symmetries"
  DESTINATION ${INCLUDE_INSTALL_DIR}
  FILES_MATCHING PATTERN "*.h"
)

# Install ...Config.cmake and ...ConfigVersion into specified directory
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/chemical_symmetries/chemical_symmetriesConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/chemical_symmetries/chemical_symmetriesConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
  COMPONENT ${ComponentName}
)

# Source groups for Visual Studio
source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
