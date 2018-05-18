# import assignment and symmetry libraries so we can specify a dependency
find_package(stereopermutation 0.1.0 REQUIRED)
find_package(chemical_symmetries 0.1.0 REQUIRED)
find_package(temple 0.1.0 REQUIRED)
find_package(cyclic_polygons 0.1.0 REQUIRED)
find_package(scine COMPONENTS Delib REQUIRED)
find_package(BLAS REQUIRED)
find_package(dlib REQUIRED)

set(Boost_USE_STATIC_LIBS OFF)
find_package(Boost REQUIRED COMPONENTS filesystem system)

# add include directories
include_directories(SYSTEM ${RDL_INCLUDE_DIR} include)

set(molassembler_VERSION 0.1.0)

# collect headers
set(HEADERS_all
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceBoundsMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceGeometry.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/Error.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ExplicitGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/Gor1.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraph.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraphBoost.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraphHelpers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MetricMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MoleculeSpatialModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementDebugData.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementProblem.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ValueBounds.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/dlibAdaptors.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/dlibDebugAdaptors.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/generateConformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AdjacencyMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AnalysisHelpers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AtomInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondDistance.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BoundsFromSymmetry.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BuildTypeSwitch.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/CNStereocenter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/CommonTrig.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/CycleData.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DelibHelpers.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/EZStereocenter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphAlgorithms.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphDistanceMatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IterateStereocenterPermutations.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/LocalGeometryModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Log.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/MemFnCache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/MolGraphWriter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/OrderDiscoveryHelper.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/PersistentIndexMoleculeInterface.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RangeForTemporary.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingTree.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/SequenceSet.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StdlibTypeAlgorithms.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereocenter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereocenterList.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Version.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/common_typedefs.h
)

# collect sources
set(SOURCES_all
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/DistanceBoundsMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ExplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ImplicitGraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MetricMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/MoleculeSpatialModel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/RefinementProblem.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/ValueBounds.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DistanceGeometry/generateConformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AdjacencyMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AtomInfo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondDistance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/CNStereocenter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/CommonTrig.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/CycleData.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/DelibHelpers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/EZStereocenter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphDistanceMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/IO.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/LocalGeometryModel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Log.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/MolGraphWriter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingTree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StdlibTypeAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereocenter.cpp
)

# Static variant
add_library(molassemblerStatic STATIC ${SOURCES_all} ${HEADERS_all})
target_include_directories(
  molassemblerStatic PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/molassembler>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/molassembler>
    ${EIGEN3_INCLUDE_DIRS}
)

# Shared variant
add_library(molassembler SHARED ${SOURCES_all} ${HEADERS_all})
target_include_directories(
  molassembler PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/molassembler>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/molassembler>
    ${EIGEN3_INCLUDE_DIRS}
)

# Compilation options
if(UNIX)
  add_definitions(-fPIC)
endif()

target_compile_options(molassemblerStatic PRIVATE -Wall -Wpedantic -Wextra)
target_compile_options(molassembler PRIVATE -Wall -Wpedantic -Wextra)

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
    scine::Delib
    stereopermutation
    boost_outcome
    chemical_symmetries
    RingDecomposerLibStatic
    gor1
    dlib::dlib
    ${BLAS_LIBRARIES}
  PRIVATE
    cyclic_polygons
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
)

target_link_libraries(
  molassemblerStatic
  PUBLIC
    Boost::boost
    scine::DelibStatic
    stereopermutation
    boost_outcome
    chemical_symmetries
    RingDecomposerLibStatic
    gor1
    dlib::dlib
    ${BLAS_LIBRARIES}
  PRIVATE
    cyclic_polygons
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
