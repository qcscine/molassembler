# Dependencies
if(MOLASSEMBLER_PARALLELIZE)
  find_package(OpenMP)
endif()

if(NOT TARGET Boost::filesystem OR NOT TARGET Boost::system)
  find_package(Boost REQUIRED COMPONENTS filesystem system)
endif()

set(_LIBPREFIX ${CMAKE_CURRENT_SOURCE_DIR}/molassembler)

set(HEADERS_all
  ${_LIBPREFIX}/AngstromPositions.h
  ${_LIBPREFIX}/AtomStereopermutator.h
  ${_LIBPREFIX}/BondOrders.h
  ${_LIBPREFIX}/BondStereopermutator.h
  ${_LIBPREFIX}/Conformers.h
  ${_LIBPREFIX}/Containers/OrderDiscoveryHelper.h
  ${_LIBPREFIX}/Cycles.h
  ${_LIBPREFIX}/Descriptors.h
  ${_LIBPREFIX}/Detail/BuildTypeSwitch.h
  ${_LIBPREFIX}/Detail/Cartesian.h
  ${_LIBPREFIX}/DirectedConformerGenerator.h
  ${_LIBPREFIX}/DistanceGeometry/ConformerGeneration.h
  ${_LIBPREFIX}/DistanceGeometry/DirectedConformerGeneratorImpl.h
  ${_LIBPREFIX}/DistanceGeometry/DistanceBoundsMatrix.h
  ${_LIBPREFIX}/DistanceGeometry/DistanceGeometry.h
  ${_LIBPREFIX}/DistanceGeometry/Error.h
  ${_LIBPREFIX}/DistanceGeometry/ExplicitBoundsGraph.h
  ${_LIBPREFIX}/DistanceGeometry/Gor1.h
  ${_LIBPREFIX}/DistanceGeometry/ImplicitBoundsGraphBoost.h
  ${_LIBPREFIX}/DistanceGeometry/ImplicitBoundsGraph.h
  ${_LIBPREFIX}/DistanceGeometry/ImplicitBoundsGraphHelpers.h
  ${_LIBPREFIX}/DistanceGeometry/MetricMatrix.h
  ${_LIBPREFIX}/DistanceGeometry/RefinementMeta.h
  ${_LIBPREFIX}/DistanceGeometry/SpatialModel.h
  ${_LIBPREFIX}/DistanceGeometry/TetrangleSmoothing.h
  ${_LIBPREFIX}/DistanceGeometry/ValueBounds.h
  ${_LIBPREFIX}/Editing.h
  ${_LIBPREFIX}/Graph.h
  ${_LIBPREFIX}/Graph/Bridge.h
  ${_LIBPREFIX}/Graph/Canonicalization.h
  ${_LIBPREFIX}/Graph/GraphAlgorithms.h
  ${_LIBPREFIX}/Graph/PrivateGraph.h
  ${_LIBPREFIX}/GraphAlgorithms.h
  ${_LIBPREFIX}/IO.h
  ${_LIBPREFIX}/IO/Base64.h
  ${_LIBPREFIX}/IO/BinaryHandler.h
  ${_LIBPREFIX}/IO/SmilesMoleculeBuilder.h
  ${_LIBPREFIX}/IO/SmilesParseData.h
  ${_LIBPREFIX}/IO/SmilesParser.h
  ${_LIBPREFIX}/Interpret.h
  ${_LIBPREFIX}/Isomers.h
  ${_LIBPREFIX}/Log.h
  ${_LIBPREFIX}/Modeling/AtomInfo.h
  ${_LIBPREFIX}/Modeling/BondDistance.h
  ${_LIBPREFIX}/Modeling/CommonTrig.h
  ${_LIBPREFIX}/Modeling/ShapeInference.h
  ${_LIBPREFIX}/Molecule.h
  ${_LIBPREFIX}/Molecule/AtomEnvironmentHash.h
  ${_LIBPREFIX}/Molecule/MolGraphWriter.h
  ${_LIBPREFIX}/Molecule/MoleculeImpl.h
  ${_LIBPREFIX}/Molecule/RankingTree.h
  ${_LIBPREFIX}/Options.h
  ${_LIBPREFIX}/Prng.h
  ${_LIBPREFIX}/RankingInformation.h
  ${_LIBPREFIX}/Serialization.h
  ${_LIBPREFIX}/StereopermutatorList.h
  ${_LIBPREFIX}/Stereopermutators/AbstractPermutations.h
  ${_LIBPREFIX}/Stereopermutators/AtomStereopermutatorImpl.h
  ${_LIBPREFIX}/Stereopermutators/BondStereopermutatorImpl.h
  ${_LIBPREFIX}/Stereopermutators/CycleFeasibility.h
  ${_LIBPREFIX}/Stereopermutators/FeasiblePermutations.h
  ${_LIBPREFIX}/Stereopermutators/ShapeVertexMaps.h
  ${_LIBPREFIX}/Subgraphs.h
  ${_LIBPREFIX}/Types.h
)

# collect sources
set(SOURCES_all
  ${_LIBPREFIX}/AngstromPositions.cpp
  ${_LIBPREFIX}/AtomStereopermutator.cpp
  ${_LIBPREFIX}/BondOrders.cpp
  ${_LIBPREFIX}/BondStereopermutator.cpp
  ${_LIBPREFIX}/Conformers.cpp
  ${_LIBPREFIX}/Cycles.cpp
  ${_LIBPREFIX}/Descriptors.cpp
  ${_LIBPREFIX}/Detail/Cartesian.cpp
  ${_LIBPREFIX}/DirectedConformerGenerator.cpp
  ${_LIBPREFIX}/DistanceGeometry/ConformerGeneration.cpp
  ${_LIBPREFIX}/DistanceGeometry/DirectedConformerGeneratorImpl.cpp
  ${_LIBPREFIX}/DistanceGeometry/DistanceBoundsMatrix.cpp
  ${_LIBPREFIX}/DistanceGeometry/ExplicitBoundsGraph.cpp
  ${_LIBPREFIX}/DistanceGeometry/ImplicitBoundsGraph.cpp
  ${_LIBPREFIX}/DistanceGeometry/MetricMatrix.cpp
  ${_LIBPREFIX}/DistanceGeometry/SpatialModel.cpp
  ${_LIBPREFIX}/DistanceGeometry/TetrangleSmoothing.cpp
  ${_LIBPREFIX}/DistanceGeometry/ValueBounds.cpp
  ${_LIBPREFIX}/Editing.cpp
  ${_LIBPREFIX}/Graph.cpp
  ${_LIBPREFIX}/Graph/Canonicalization.cpp
  ${_LIBPREFIX}/Graph/GraphAlgorithms.cpp
  ${_LIBPREFIX}/Graph/PrivateGraph.cpp
  ${_LIBPREFIX}/GraphAlgorithms.cpp
  ${_LIBPREFIX}/GraphIterators.cpp
  ${_LIBPREFIX}/IO.cpp
  ${_LIBPREFIX}/IO/Base64.cpp
  ${_LIBPREFIX}/IO/BinaryHandler.cpp
  ${_LIBPREFIX}/IO/SmilesMoleculeBuilder.cpp
  ${_LIBPREFIX}/IO/SmilesParser.cpp
  ${_LIBPREFIX}/Interpret.cpp
  ${_LIBPREFIX}/Isomers.cpp
  ${_LIBPREFIX}/Log.cpp
  ${_LIBPREFIX}/Modeling/AtomInfo.cpp
  ${_LIBPREFIX}/Modeling/BondDistance.cpp
  ${_LIBPREFIX}/Modeling/CommonTrig.cpp
  ${_LIBPREFIX}/Modeling/ShapeInference.cpp
  ${_LIBPREFIX}/Molecule.cpp
  ${_LIBPREFIX}/Molecule/AtomEnvironmentHash.cpp
  ${_LIBPREFIX}/Molecule/MolGraphWriter.cpp
  ${_LIBPREFIX}/Molecule/MoleculeImpl.cpp
  ${_LIBPREFIX}/Molecule/RankingTree.cpp
  ${_LIBPREFIX}/Options.cpp
  ${_LIBPREFIX}/Prng.cpp
  ${_LIBPREFIX}/RankingInformation.cpp
  ${_LIBPREFIX}/Serialization.cpp
  ${_LIBPREFIX}/StereopermutatorList.cpp
  ${_LIBPREFIX}/Stereopermutators/AbstractPermutations.cpp
  ${_LIBPREFIX}/Stereopermutators/AtomStereopermutatorImpl.cpp
  ${_LIBPREFIX}/Stereopermutators/BondStereopermutatorImpl.cpp
  ${_LIBPREFIX}/Stereopermutators/CycleFeasibility.cpp
  ${_LIBPREFIX}/Stereopermutators/FeasiblePermutations.cpp
  ${_LIBPREFIX}/Stereopermutators/ShapeVertexMaps.cpp
  ${_LIBPREFIX}/Subgraphs.cpp
  ${_LIBPREFIX}/Types.cpp
)

unset(_LIBPREFIX)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
  @ONLY
)

# Object library to avoid compiling twice for shared and static variants
add_library(molassembler_obj OBJECT
  ${SOURCES_all}
  ${HEADERS_all}
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
)

include(GenerateExportHeader)
generate_export_header(molassembler_obj
  BASE_NAME masm
  EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Export.h
)

set_target_properties(molassembler_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(molassembler_obj PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

# Object library build requirements
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  # this is clearly the preferred method of handling this mess, but only
  # available from 3.12 onward. The proper linking an interface population is
  # handled below
  target_link_libraries(molassembler_obj PRIVATE
    temple
    Boost::boost
    ${MOLASSEMBLER_UTILS_TARGET}
    boost_outcome
    RingDecomposerLibStatic
    stereopermutation
    cyclic_polygons
    gor1
    json
    Boost::filesystem
    Boost::system
    nauty
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )
else()
  target_include_directories(molassembler_obj
    PUBLIC
    $<TARGET_PROPERTY:${MOLASSEMBLER_UTILS_TARGET},INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:boost_outcome,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:RingDecomposerLibStatic,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:Boost::boost,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:shapes,INTERFACE_INCLUDE_DIRECTORIES>
    PRIVATE
      $<TARGET_PROPERTY:temple,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:stereopermutation,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:cyclic_polygons,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:gor1,INTERFACE_INCLUDE_DIRECTORIES>
  )
endif()
target_link_eigen(molassembler_obj PRIVATE)
target_compile_options(molassembler_obj PRIVATE
  ${MOLASSEMBLER_CXX_FLAGS}
  $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
)
# Suppress warnings coming from various external libraries
target_include_directories(molassembler_obj
  SYSTEM PRIVATE
    ${Boost_INCLUDE_DIR}
    ${RDL_INCLUDE_DIR}
    $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:json,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:nauty,INTERFACE_INCLUDE_DIRECTORIES>
)
set_target_properties(molassembler_obj PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)

# Shared library
add_library(molassembler SHARED
  $<TARGET_OBJECTS:molassembler_obj>
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Export.h
)
target_include_directories(molassembler
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

add_library(molassembler::molassembler ALIAS molassembler)

# Linking for the shared version of molassembler is fairly straightforward.
# Some notes:
# - Boost::boost Types like boost::optional form part of the public interface
# - Boost::filesystem, ::system are required for implementation details and
#   hence PRIVATE, subsumed since we are using the boost static libraries
# - boost_outcome Result types form part of the public interface, distributed
#   as a COMPONENT in this build tree (see extern)
# - All other PRIVATE libraries are STATIC or INTERFACE libraries that we can
#   subsume
target_link_libraries(molassembler
  PUBLIC
    Boost::boost
    ${MOLASSEMBLER_UTILS_TARGET}
    boost_outcome
    shapes
  PRIVATE
    temple
    stereopermutation
    cyclic_polygons
    gor1
    json
    Boost::filesystem
    Boost::system
    RingDecomposerLibStatic
    nauty
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)
target_link_eigen(molassembler PUBLIC)
set_target_properties(molassembler PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)

# Static library
add_library(molassembler_static STATIC
  $<TARGET_OBJECTS:molassembler_obj>
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Export.h
)
set_target_properties(molassembler_static PROPERTIES
  OUTPUT_NAME "molassembler"
)
target_include_directories(molassembler_static PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
add_library(molassembler::molassembler_static ALIAS molassembler_static)

# The linking for the static version of molassembler is more complicated.
#
# Regarding the differentiation between public and private dependencies, it is
# important to note two things:
# - When B is linked into A, and both A and B are static libraries, no
#   transfer of symbols from B into A is performed at all. If you wish to
#   perform something of that likeness, you need to literally concatenate the
#   two static libraries together with the platform archiver. This does NOT
#   perform anything of the normal duties of a linker in the sense that no
#   distinction is made whether the symbols in B are actually needed or not,
#   the archiver merely pastes these libraries together.
# - For the above reason, and since static libraries carry no information by
#   themselves about where encountered undefined symbols may be found, CMake
#   handles PRIVATE linking of a static library as if it were PUBLIC and
#   propagates the linking requirements of B onto A.
#
# In effect, no distinction is made between PUBLIC and PRIVATE for static A
# and B. The behavioral claims here are based on
# https://cmake.org/pipermail/cmake/2016-May/063400.html
#
# However, some of the dependencies of this static version of molassembler are
# targets declared via add_library(x INTERFACE), i.e. they are header-only. We
# need them only for implementation details and they do not provide any
# symbols, so they are PRIVATE. CMake still references their targets as
# LINK_ONLY dependencies in INTERFACE_LINK_LIBRARIES for the static library!
#
# Yet we want only the include information from these libraries so that
# molassembler_static can be built properly, yet not the link dependency.
#
# So we link the special cases, populating all dependency information:
target_link_libraries(molassembler_static PRIVATE
  cyclic_polygons
  gor1
  json
)

# And then reset the libraries listed in INTERFACE_LINK_LIBRARIES (which are
# only the aforementioned libraries specified in a LINK_ONLY generator
# expression). This removes the link dependency but does retain the
# information required to build molassembler_static from those targets.
set_target_properties(molassembler_static PROPERTIES
  INTERFACE_LINK_LIBRARIES ""
)

# And finally we link in the rest of the public interface libraries
target_link_libraries(molassembler_static PUBLIC
  Boost::boost
  ${MOLASSEMBLER_UTILS_TARGET}
  boost_outcome
  RingDecomposerLibStatic
  stereopermutation
  Boost::filesystem
  Boost::system
  nauty
  $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)
target_link_eigen(molassembler_static PUBLIC)
set_target_properties(molassembler_static PROPERTIES
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)

install(
  TARGETS molassembler molassembler_static
  EXPORT molassemblerTargets
  DESTINATION lib
)

molassembler_install_component_cmake_files(
  COMPONENT molassembler
  EXPORT_NAME molassemblerTargets
)

# We wish to install headers only and not install any subdirectories since
# molassembler's public API is defined via the top-level-directory headers
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/molassembler"
  DESTINATION ${INCLUDE_INSTALL_DIR}
  COMPONENT molassembler
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*.hxx"
)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
    ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Export.h
  DESTINATION ${INCLUDE_INSTALL_DIR}/molassembler
  COMPONENT molassembler
)

source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
