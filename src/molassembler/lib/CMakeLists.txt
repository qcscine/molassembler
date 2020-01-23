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
  ${_LIBPREFIX}/DistanceGeometry/DistanceGeometry.cpp
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
    Scine::UtilsOS
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
      $<TARGET_PROPERTY:Scine::UtilsOS,INTERFACE_INCLUDE_DIRECTORIES>
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
target_include_directories(molassembler_obj SYSTEM PRIVATE
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

function(molassembler_includes_and_properties target_name)
  target_include_directories(${target_name} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  )
  set_target_properties(${target_name} PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
    OUTPUT_NAME "molassembler"
  )
endfunction()

function(molassembler_library_links target_name)
  # Some of the dependencies of this possibly static version of molassembler
  # are INTERFACE libraries. We need them only for implementation details and
  # they do not provide any symbols, so they are linked PRIVATE. CMake still
  # references their targets as LINK_ONLY dependencies in
  # INTERFACE_LINK_LIBRARIES for the static library!
  #
  # Yet we want only the include information from these libraries so that
  # the library can be built properly, yet not the link dependency.
  # See in particular: https://gitlab.kitware.com/cmake/cmake/issues/15415
  target_include_directories(${target_name} PRIVATE
    $<TARGET_PROPERTY:cyclic_polygons,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:gor1,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:json,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:temple,INTERFACE_INCLUDE_DIRECTORIES>
  )

  # And then we link in the rest of the libraries just like in the shared case
  target_link_libraries(${target_name}
    PUBLIC
      Scine::UtilsOS
      boost_outcome
      shapes
    PRIVATE
      stereopermutation
      Boost::filesystem
      Boost::system
      RingDecomposerLibStatic
      nauty
      $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )

  target_link_eigen(${target_name} PUBLIC)
endfunction()

# Main library, of type as determined by BUILD_SHARED_LIBS
add_library(molassembler
  $<TARGET_OBJECTS:molassembler_obj>
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Export.h
)
molassembler_includes_and_properties(molassembler)
molassembler_library_links(molassembler)
add_library(molassembler::molassembler ALIAS molassembler)

if(BUILD_SHARED_LIBS)
  # If molassembler is shared, we still need a static variant for linking into
  # the tests and other binaries that hook deeper into molassembler than the
  # public API
  add_library(molassembler_static STATIC
    $<TARGET_OBJECTS:molassembler_obj>
    ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
    ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Export.h
  )
  molassembler_includes_and_properties(molassembler_static)
  molassembler_library_links(molassembler_static)
  add_library(molassembler::molassembler_static ALIAS molassembler_static)

  # Target name for others to hook into where needed
  set(MOLASSEMBLER_STATIC_TARGET molassembler_static PARENT_SCOPE)
else()
  # If molassembler is static, then that's all we need
  set(MOLASSEMBLER_STATIC_TARGET molassembler PARENT_SCOPE)
endif()

install(
  TARGETS molassembler
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
