# Dependencies
if(MOLASSEMBLER_PARALLELIZE)
  find_package(OpenMP)
endif()

if(NOT TARGET Boost::filesystem OR NOT TARGET Boost::system)
  find_package(Boost REQUIRED COMPONENTS filesystem system)
endif()

# Add RingDecomposerLib include directory
include_directories(SYSTEM ${RDL_INCLUDE_DIR} include)

set(HEADERS_all
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AngstromWrapper.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AtomStereopermutator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondOrders.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondStereopermutator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Conformers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Containers/OrderDiscoveryHelper.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Cycles.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Descriptors.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/BuildTypeSwitch.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/Cartesian.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DirectedConformerGenerator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ConformerGeneration.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DirectedConformerGeneratorImpl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceBoundsMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceGeometry.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/Error.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ExplicitGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/Gor1.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraphBoost.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraphHelpers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MetricMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementDebugData.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementMeta.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/SpatialModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/TetrangleSmoothing.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ValueBounds.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Editing.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/Bridge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/Canonicalization.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/GraphAlgorithms.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/InnerGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphAlgorithms.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Interpret.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Isomers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO/Base64.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO/BinaryHandler.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Log.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/AtomInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/BondDistance.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/CommonTrig.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/LocalGeometryModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/AtomEnvironmentHash.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/MolGraphWriter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/MoleculeImpl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/RankingTree.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Options.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/OuterGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Patterns.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/PRNG.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Serialization.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereopermutatorList.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/AbstractPermutations.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/AtomStereopermutatorImpl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/BondStereopermutatorImpl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/CycleFeasibility.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/FeasiblePermutations.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/ShapeVertexMaps.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Subgraphs.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Types.h
)

# collect sources
set(SOURCES_all
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AngstromWrapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AtomStereopermutator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondOrders.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondStereopermutator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Conformers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Cycles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Descriptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/Cartesian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DirectedConformerGenerator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ConformerGeneration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DirectedConformerGeneratorImpl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceBoundsMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ExplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MetricMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/SpatialModel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/TetrangleSmoothing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ValueBounds.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Editing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/Canonicalization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/GraphAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/InnerGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Interpret.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Isomers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO/BinaryHandler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Log.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/AtomInfo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/BondDistance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/CommonTrig.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/LocalGeometryModel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/AtomEnvironmentHash.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/MolGraphWriter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/MoleculeImpl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/RankingTree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Options.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/OuterGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/OuterGraphIterators.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Patterns.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/PRNG.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Serialization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereopermutatorList.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/AbstractPermutations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/AtomStereopermutatorImpl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/BondStereopermutatorImpl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/CycleFeasibility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/FeasiblePermutations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/ShapeVertexMaps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Subgraphs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Types.cpp
)

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
  @ONLY
)

# Shared library
add_library(molassembler SHARED
  ${SOURCES_all}
  ${HEADERS_all}
  ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
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
target_link_libraries(
  molassembler
  PUBLIC
    Boost::boost
    Scine::UtilsOS
    boost_outcome
    RingDecomposerLibStatic
    shapes
  PRIVATE
    temple
    stereopermutation
    cyclic_polygons
    gor1
    json
    spectra
    Boost::filesystem
    Boost::system
    nautyStatic
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)
add_eigen(molassembler)

# Compiler flags and definitions
target_compile_options(molassembler PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
target_compile_definitions(molassembler PRIVATE BOOST_ALL_DYN_LINK)

# Suppress warnings coming from various external libraries
target_include_directories(molassembler SYSTEM PRIVATE
  ${Boost_INCLUDE_DIR}
  $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:json,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:spectra,INTERFACE_INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:nautyStatic,INTERFACE_INCLUDE_DIRECTORIES>
)

# -- END Linking

install(
  TARGETS molassembler
  EXPORT molassemblerTargets
  DESTINATION $<$<NOT:$<CONFIG:Release>>:$<CONFIG>/>lib
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
  FILES ${CMAKE_CURRENT_BINARY_DIR}/molassembler/Version.h
  DESTINATION ${INCLUDE_INSTALL_DIR}
  COMPONENT molassembler
)

source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
