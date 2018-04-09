# TODO must get rid of absolute paths here!
set(RDL_LIBRARY_DIR "/usr/local/lib/")
set(RDL_INCLUDE_DIR "/usr/local/include/")

# dependencies: libraries
find_library(RDL_LIBRARY RingDecomposerLibStatic ${RDL_LIBRARY_DIR})

# import assignment and symmetry libraries so we can specify a dependency
find_package(stereopermutation 0.1.0 REQUIRED)
find_package(chemical_symmetries 0.1.0 REQUIRED)
find_package(temple 0.1.0 REQUIRED)
find_package(cyclic_polygons 0.1.0 REQUIRED)

add_library(rdl STATIC IMPORTED)
set_property(TARGET rdl PROPERTY IMPORTED_LOCATION ${RDL_LIBRARY})

# add include directories
include_directories(SYSTEM ${RDL_INCLUDE_DIR} include)

set(MOLASSEMBLER_VERSION 0.1.0)

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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Edges.h
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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/LocalGeometryModel.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Log.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/MolGraphWriter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingTree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StdlibTypeAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereocenter.cpp
)

# specify library target
add_library(molassembler STATIC ${SOURCES_all} ${HEADERS_all})
target_include_directories(
  molassembler
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/molassembler>
)
target_include_directories(
  molassembler
    SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/molassembler>
)

# build includes
target_include_directories(molassembler PUBLIC ${EIGEN3_INCLUDE_DIRS})

# Compilation options
if(UNIX)
  add_definitions(-fPIC)
endif()

target_compile_options(molassembler PRIVATE -Wall -Wpedantic -Wextra)

# Switch debug and release build types
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  message(STATUS "Selected release build")
  # no debug information, disable assert
  target_compile_definitions(molassembler PRIVATE NDEBUG)
  target_compile_options(molassembler PUBLIC -O3 -march=native)
else() # assume Debug
  message(STATUS "Assume debug build")
  target_compile_options(molassembler PUBLIC -g)
  # target_compile_options(molassembler PUBLIC -g -fsanitize=address -fsanitize=undefined)
  # set_target_properties(molassembler PROPERTIES LINK_FLAGS "-fsanitize=address -fsanitize=undefined")
endif()

# enable clang-tidy warnings
#if(CLANG_TIDY_EXE)
#  set_target_properties(
#    molassembler PROPERTIES
#    CXX_CLANG_TIDY "${DO_CLANG_TIDY}"
#  )
#endif()

# Linking
target_link_libraries(
  molassembler
  PUBLIC
    Delib::DelibStatic
    stereopermutation
    cyclic_polygons
    outcome
    gor1
    chemicalsymmetries
    blas
    rdl
)

# Installation
set(ConfigPackageLocation lib/cmake/molassembler)
set(IncludeInstallDir include/)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/molassemblerConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/molassembler/molassemblerConfig.cmake"
  INSTALL_DESTINATION ${ConfigPackageLocation}
  PATH_VARS IncludeInstallDir
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/molassembler/molassemblerConfigVersion.cmake"
  VERSION ${MOLASSEMBLER_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(
  TARGETS molassembler
  EXPORT molassemblerTargets
  DESTINATION lib
)

install(
  EXPORT molassemblerTargets
  DESTINATION "${ConfigPackageLocation}"
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/molassembler"
  DESTINATION IncludeInstallDir
  FILES_MATCHING
  PATTERN "*.h"
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/molassembler/molassemblerConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/molassembler/molassemblerConfigVersion.cmake"
  DESTINATION ${ConfigPackageLocation}
)

source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
