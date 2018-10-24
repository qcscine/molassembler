# declaration
project(GORProject)

set(gor1_VERSION 0.1.0)

# dependencies
if(NOT TARGET Boost::unit_test_framework)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework)
endif()

# Files
set(HEADERS_All ${CMAKE_CURRENT_SOURCE_DIR}/gor1/Gor1.h)

enable_testing()

# Interface library
add_library(gor1 INTERFACE)
target_sources(
  gor1 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/gor1/Gor1.h>
)
target_include_directories(
  gor1 INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)
target_include_directories(
  gor1 SYSTEM INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

target_compile_options(gor1 INTERFACE ${MOLASSEMBLER_CXX_FLAGS})

# Installation
install(
  TARGETS gor1
  EXPORT gor1Targets
  DESTINATION ${INCLUDE_INSTALL_DIR}
)

molassembler_install_component_cmake_files(
  COMPONENT gor1
  EXPORT_NAME gor1Targets
)

molassembler_install_component_headers(
  COMPONENT gor1
  DIRECTORY gor1
)

# Tests
add_executable(gor1Tests test.cpp)
target_link_libraries(gor1Tests PRIVATE Boost::boost Boost::unit_test_framework)

add_test(NAME gor1Tests COMMAND gor1Tests)

set_property(TARGET gor1Tests PROPERTY FOLDER "Tests")

# Source groups
source_group("Headers" FILES ${HEADERS_All})
