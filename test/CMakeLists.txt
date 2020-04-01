if(NOT SCINE_BUILD_TESTS)
  return()
endif()

if(NOT TARGET Boost::unit_test_framework
    OR NOT TARGET Boost::regex
    OR NOT TARGET Boost::filesystem
)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework regex filesystem QUIET)
endif()

set(MainTestFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/AtomStereopermutator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/BondStereopermutator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cartesian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Coherence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/CommonTrig.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cycles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Descriptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Editing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Fixtures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/PrivateGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/OrderDiscoveryHelper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingSpotCheck.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Serialization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Smiles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Subgraphs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Vsepr.cpp
)

set(DGTestFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/DirectedConformerGenerator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/ExplicitBoundsGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/FixedPositions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/GraphShortestDistances.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/ImplicitBoundsGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/MetricMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/RefinementProblem.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/Smoothing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DistanceGeometry/main.cpp
)

# Testing targets
add_executable(MolassemblerTests ${MainTestFiles})
target_link_libraries(MolassemblerTests
  PRIVATE
    Boost::boost
    Boost::regex
    gor1
    ${MOLASSEMBLER_STATIC_TARGET}
    stereopermutation
    Boost::unit_test_framework
    Boost::filesystem
)
target_compile_options(MolassemblerTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_include_directories(MolassemblerTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME MASMMainTests
  COMMAND MolassemblerTests
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../data/"
)

add_executable(DistanceGeometryTests ${DGTestFiles})
target_link_libraries(DistanceGeometryTests
  PRIVATE
    Boost::boost
    Boost::regex
    gor1
    ${MOLASSEMBLER_STATIC_TARGET}
    stereopermutation
    Boost::unit_test_framework
    Boost::filesystem
)
include(LibType)
determine_lib_type(Boost::unit_test_framework)
if(${LIB_TYPE} STREQUAL "SHARED_LIBRARY")
  target_compile_definitions(MolassemblerTests PRIVATE BOOST_ALL_DYN_LINK)
  target_compile_definitions(DistanceGeometryTests PRIVATE BOOST_ALL_DYN_LINK)
endif()

target_compile_options(DistanceGeometryTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_include_directories(DistanceGeometryTests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
add_test(
  NAME MASMDistanceGeometryTests
  COMMAND DistanceGeometryTests
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../data/"
)
set_property(TARGET DistanceGeometryTests PROPERTY FOLDER "Tests")
