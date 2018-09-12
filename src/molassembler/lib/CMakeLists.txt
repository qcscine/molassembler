# import assignment and symmetry libraries so we can specify a dependency
find_package(stereopermutation 0.1.0 REQUIRED)
find_package(chemical_symmetries 0.1.0 REQUIRED)
find_package(temple 0.1.0 REQUIRED)
find_package(cyclic_polygons 0.1.0 REQUIRED)
find_package(Delib REQUIRED)
find_package(BLAS REQUIRED)
find_package(dlib REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system)

# add include directories
include_directories(SYSTEM ${RDL_INCLUDE_DIR} include)

set(molassembler_VERSION 0.1.0)

# collect headers
# if a big update happens, use $ find . -name "*.h" | sort
set(HEADERS_all
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AngstromWrapper.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AtomStereocenter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondOrders.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondStereocenter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Conformers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Containers/OrderDiscoveryHelper.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Cycles.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Descriptors.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/AnalysisHelpers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/BuildTypeSwitch.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/DelibHelpers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/StdlibTypeAlgorithms.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ConformerGeneration.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceBoundsMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceGeometry.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/dlibAdaptors.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/dlibDebugAdaptors.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/Error.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ExplicitGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/Gor1.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraphBoost.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraphHelpers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MetricMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementDebugData.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementProblem.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/SpatialModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ValueBounds.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/Bridge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/GraphAlgorithms.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/InnerGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Interpret.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO/Base64.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO/FileHandlers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Log.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/AtomInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/BondDistance.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/CommonTrig.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/LocalGeometryModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/AtomEnvironmentHash.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/MolGraphWriter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/RankingTree.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Options.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/OuterGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Serialization.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereocenterList.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereocenters/PermutationState.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Version.h
)

# collect sources
set(SOURCES_all
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AngstromWrapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AtomStereocenter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondOrders.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondStereocenter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Conformers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Cycles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Descriptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/DelibHelpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Detail/StdlibTypeAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ConformerGeneration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceBoundsMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ExplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MetricMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementProblem.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/SpatialModel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ValueBounds.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/GraphAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Graph/InnerGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Interpret.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO/FileHandlers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Log.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/AtomInfo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/BondDistance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/CommonTrig.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Modeling/LocalGeometryModel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/AtomEnvironmentHash.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/MolGraphWriter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule/RankingTree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Options.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/OuterGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/OuterGraphIterators.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Serialization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereocenterList.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereocenters/PermutationState.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Types.cpp
)

# Static variant
add_library(molassemblerStatic STATIC ${SOURCES_all} ${HEADERS_all})
target_include_directories(molassemblerStatic
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
target_include_directories(molassemblerStatic
  SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIRS}
)

# Shared variant
add_library(molassembler SHARED ${SOURCES_all} ${HEADERS_all})
target_include_directories(molassembler
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
target_include_directories(molassembler
  SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIRS}
)

# Linking
# NOTE: For the shared variant, it is necessary to use the cmake variables
# instead of the imported targets since the imported targets are not changed
# between find_package(Boost) invocations, but the Boost scoped variables are.
# Otherwise, Boost::filesystem is the static library (despite setting
# Boost_USE_STATIC_LIBS OFF) above.
target_link_libraries(
  molassembler
  PUBLIC
    Boost::boost
    delib::Delib
    boost_outcome
    chemical_symmetries
    RingDecomposerLibStatic
    dlib::dlib
    ${BLAS_LIBRARIES}
  PRIVATE
    cyclic_polygons
    stereopermutation
    gor1
    json
    Boost::filesystem
    Boost::system
)

target_link_libraries(
  molassemblerStatic
  PUBLIC
    Boost::boost
    delib::DelibStatic
    boost_outcome
    chemical_symmetries
    RingDecomposerLibStatic
    dlib::dlib
    ${BLAS_LIBRARIES}
  PRIVATE
    cyclic_polygons
    stereopermutation
    gor1
    json
    Boost::filesystem
    Boost::system
)

install(
  TARGETS molassembler molassemblerStatic
  EXPORT molassemblerTargets
  DESTINATION lib
)

molassembler_install_component_cmake_files(
  COMPONENT molassembler
  EXPORT_NAME molassemblerTargets
)

molassembler_install_component_headers(
  COMPONENT molassembler
  DIRECTORY molassembler
)

source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
