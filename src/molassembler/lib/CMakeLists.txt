# Dependencies
find_package(dlib REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem system)

# For CI purposes, it is very useful to just download in-tree dependencies from
# GitLab. This avoids having to install the dependencies in the docker image and
# also avoids recompiling the dependencies in every build due to the CI cache.
# However, downloading the dependencies adds the files into the tree, but does
# NOT mean that the dependencies are also installed!
#
# Every human user of this project should independently install each project
# and specify delib_DIR or CMAKE_PREFIX_PATH
if(MOLASSEMBLER_DOWNLOAD_DEPENDENCIES)
  # Attempt to find Delib locally
  find_package(Delib QUIET)
  if(NOT TARGET delib::Delib)
    message(STATUS "Could not find Delib locally, downloading from gitlab.")
    include(DelibImport)
    import_delib()
    if(NOT TARGET delib::Delib)
      message(FATAL_ERROR "Delib could not be downloaded from gitlab!")
    endif()
  endif()
else()
  find_package(Delib REQUIRED)
endif()

# Some notes on BLAS / LAPACK:
# Both dlib and Eigen can use BLAS / LAPACK routines.
# - Eigen, as a header-only library, merely requires use to add the definitions
#   and link against the library.
# - dlib sets DLIB_USE_BLAS and DLIB_USE_LAPACK for itself at compile-time if
#   it can find either type of library in two fashions:
#   - cblas bindings are undefined in the static library, but the dlib::dlib
#     target has the found libraries in its INTERFACE_LINK_LIBRARIES
#   - it configures config.h.in to reflect the use of DLIB_USE_BLAS and
#     DLIB_USE_LAPACK so that all dependent uses have the same definitions
#
# Although we enable the use of both types of libraries, currently only BLAS is
# used by either. We keep LAPACK use enabled in order to quickly add it as a
# dependency if either library can make use of it in any way.
#
# Try to find suitable BLAS and LAPACK libraries
find_package(BLAS QUIET)
# Currently, molassembler does not use LAPACK routines, so this is commented out
# find_package(LAPACK QUIET)

# Add RingDecomposerLib include directory
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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphAlgorithms.h
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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/PRNG.h
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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/GraphAlgorithms.cpp
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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/PRNG.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Serialization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereocenterList.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereocenters/PermutationState.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Types.cpp
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
target_compile_options(molassembler PRIVATE ${MOLASSEMBLER_CXX_FLAGS})

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

target_compile_options(molassemblerStatic PRIVATE ${MOLASSEMBLER_CXX_FLAGS})

# -- BEGIN Linking

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
    delib::Delib
    boost_outcome
    chemical_symmetries
    RingDecomposerLibStatic
    dlib::dlib
  PRIVATE
    stereopermutation
    cyclic_polygons
    gor1
    json
    spectra
    Boost::filesystem
    Boost::system
)

# The linking for the static version of molassembler is more complicated.
#
# Regarding the differentiation between public and private dependencies, it is
# important to note two things:
# - When B is linked into A, and both A and B are static libraries, no transfer
#   of symbols from B into A is performed at all. If you wish to perform
#   something of that likeness, you need to literally concatenate the two
#   static libraries together with the platform archiver. This does NOT perform
#   anything of the normal duties of a linker in the sense that no distinction
#   is made whether the symbols in B are actually needed or not, the archiver
#   merely pastes these libraries together.
# - For the above reason, and since static libraries carry no information by
#   themselves about where encountered undefined symbols may be found, CMake
#   handles PRIVATE linking of a static library as if it were PUBLIC and
#   propagates the linking requirements of B onto A.
#
# In effect, no distinction is made between PUBLIC and PRIVATE for static A and
# B. The behavioral claims here are based on
# https://cmake.org/pipermail/cmake/2016-May/063400.html
#
# However, some of the dependencies of this static version of molassembler are
# targets declared via add_library(x INTERFACE), i.e. they are header-only. We
# need them only for implementation details and they do not provide any symbols,
# so they are PRIVATE. CMake still references their targets as LINK_ONLY
# dependencies in INTERFACE_LINK_LIBRARIES for molassemblerStatic!
#
# Yet we want only the include information from these libraries so that
# molassemblerStatic can be built properly, yet not the link dependency.
#
# So we link the special cases, populating all dependency information:
target_link_libraries(
  molassemblerStatic
  PRIVATE
    cyclic_polygons
    gor1
    json
    spectra
)

# And then reset the libraries listed in INTERFACE_LINK_LIBRARIES (which are
# only the aforementioned libraries specified in a LINK_ONLY generator
# expression). This removes the link dependency but does retain the information
# required to build molassemblerStatic from those targets.
set_target_properties(molassemblerStatic PROPERTIES INTERFACE_LINK_LIBRARIES "")

# And finally we link in the rest of the public interface libraries
target_link_libraries(
  molassemblerStatic
  PUBLIC
    Boost::boost
    delib::DelibStatic
    boost_outcome
    chemical_symmetries
    RingDecomposerLibStatic
    stereopermutation
    dlib::dlib
    Boost::filesystem
    Boost::system
)

# Handle BLAS / LAPACK if found
if(BLAS_FOUND)
  # Add BLAS / LAPACK compile definitions
  target_compile_definitions(molassembler PUBLIC
    EIGEN_USE_BLAS
    EIGEN_USE_LAPACK
  )
  target_compile_definitions(molassemblerStatic PUBLIC
    EIGEN_USE_BLAS
    EIGEN_USE_LAPACK
  )

  target_link_libraries(molassembler PUBLIC ${BLAS_LIBARIES})
  target_link_libraries(molassemblerStatic PUBLIC ${BLAS_LIBRARIES})
endif()

# -- END Linking

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
