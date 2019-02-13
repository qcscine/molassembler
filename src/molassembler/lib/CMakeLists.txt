# Dependencies
find_package(dlib REQUIRED)

if(MOLASSEMBLER_PARALLELIZE)
  find_package(OpenMP)
endif()

if(NOT TARGET Boost::filesystem OR NOT TARGET Boost::system)
  find_package(Boost REQUIRED COMPONENTS filesystem system)
endif()

include(ScineUtils)
include(ImportUtilsOS)
import_utils_os()

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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/AtomStereopermutator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondOrders.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/BondStereopermutator.h
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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/PRNG.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Serialization.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereopermutatorList.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/AtomStereopermutatorImpl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/BondStereopermutatorImpl.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/PermutationState.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Version.h
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
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/PRNG.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/RankingInformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Serialization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/StereopermutatorList.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/AtomStereopermutatorImpl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/BondStereopermutatorImpl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Stereopermutators/PermutationState.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/molassembler/Types.cpp
)

# Object library to avoid compiling twice
add_library(molassembler_obj OBJECT ${SOURCES_all} ${HEADERS_all})
set_target_properties(molassembler_obj PROPERTIES
  POSITION_INDEPENDENT_CODE ON
)
target_include_directories(molassembler_obj PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)
# transitively populate include directories with Core's dependencies
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.12.0")
  # this is clearly the preferred method of handling this mess, but only
  # available from 3.12 onward. The proper linking an interface population is
  # handled below
  target_link_libraries(molassembler_obj
    PRIVATE
      temple
      Boost::boost
      Scine::UtilsOS
      boost_outcome
      RingDecomposerLibStatic
      dlib::dlib
      stereopermutation
      cyclic_polygons
      gor1
      json
      spectra
      Boost::filesystem
      Boost::system
      nautyStatic
  )
else()
  target_include_directories(molassembler_obj
    PUBLIC
      $<TARGET_PROPERTY:Scine::UtilsOS,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:boost_outcome,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:RingDecomposerLibStatic,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:dlib::dlib,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:Boost::boost,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:chemical_symmetries,INTERFACE_INCLUDE_DIRECTORIES>
    PRIVATE
      $<TARGET_PROPERTY:temple,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:stereopermutation,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:cyclic_polygons,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:gor1,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:json,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:spectra,INTERFACE_INCLUDE_DIRECTORIES>
      $<TARGET_PROPERTY:nautyStatic,INTERFACE_INCLUDE_DIRECTORIES>
  )
endif()
target_compile_options(molassembler_obj PRIVATE
  ${MOLASSEMBLER_CXX_FLAGS}
  $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
)
target_link_libraries(molassembler_obj PRIVATE
)
# Suppress warnings coming from various external libraries
target_include_directories(molassembler_obj
  SYSTEM PRIVATE
    ${EIGEN3_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    $<TARGET_PROPERTY:spectra,INTERFACE_INCLUDE_DIRECTORIES>
    $<TARGET_PROPERTY:nautyStatic,INTERFACE_INCLUDE_DIRECTORIES>
)

# Shared variant
add_library(molassembler SHARED $<TARGET_OBJECTS:molassembler_obj>)
target_include_directories(molassembler
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

add_library(molassembler::molassembler ALIAS molassembler)

# Static variant
add_library(molassemblerStatic STATIC $<TARGET_OBJECTS:molassembler_obj>)
target_include_directories(molassemblerStatic
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

add_library(molassembler::molassemblerStatic ALIAS molassemblerStatic)

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
    Scine::UtilsOS
    boost_outcome
    RingDecomposerLibStatic
    dlib::dlib
    chemical_symmetries
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
    Scine::UtilsOS
    boost_outcome
    RingDecomposerLibStatic
    stereopermutation
    dlib::dlib
    Boost::filesystem
    Boost::system
    nautyStatic
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
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

source_group("Headers" FILES ${HEADERS_All})
source_group("Sources" FILES ${SOURCES_All})
