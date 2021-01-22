if(NOT TARGET Boost::unit_test_framework
    OR NOT TARGET Boost::regex
    OR NOT TARGET Boost::filesystem
)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework regex filesystem QUIET)
endif()

# Copy data files to build directory for testing
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

file(GLOB_RECURSE TestSources ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
# Testing targets
add_executable(MolassemblerTests ${TestSources})
target_link_libraries(MolassemblerTests
  PRIVATE
    Boost::boost
    Boost::regex
    ${MOLASSEMBLER_STATIC_TARGET}
    Boost::unit_test_framework
    Boost::filesystem
)
target_compile_options(MolassemblerTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_include_directories(MolassemblerTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME Molassembler
  COMMAND MolassemblerTests --build_info --log_level=unit_scope
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data
)

include(LibType)
determine_lib_type(Boost::unit_test_framework)
if("${LIB_TYPE}" STREQUAL "SHARED_LIBRARY")
  target_compile_definitions(MolassemblerTests PRIVATE BOOST_ALL_DYN_LINK)
endif()
