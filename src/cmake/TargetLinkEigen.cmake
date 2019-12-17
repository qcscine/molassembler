#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
# See LICENSE.txt for details.
#

# If the target already exists, do nothing
if(NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 REQUIRED)
endif()

option(SCINE_USE_INTEL_MKL "Use the Intel MKL libraries with Eigen" ON)
option(SCINE_USE_LAPACK "Use a LAPACK library with Eigen" ON)
option(SCINE_USE_BLAS "Use a BLAS library with Eigen" ON)

# Attempt to find external linalg libraries that accelerate basic calls in
# the following order:
#
# 1. Intel MKL
# 2. LAPACK (brings BLAS too)
# 3. BLAS
#
if(NOT ADD_EIGEN_SEARCHED_EXTERNAL_LINALG_LIBRARIES)
  if(SCINE_USE_INTEL_MKL)
    include(FindMKL)
  endif()

  if(MKL_FOUND)
    find_package(OpenMP REQUIRED)
    message(STATUS "Found MKL for use with Eigen3")
  else()
    if(SCINE_USE_LAPACK)
      include(FindLAPACK)
      find_package(LAPACK QUIET)
    endif()

    if(LAPACK_FOUND)
      message(STATUS "Found LAPACK/BLAS for use with Eigen3")
    else()
      if(SCINE_USE_BLAS)
        include(FindBLAS)
        find_package(BLAS QUIET)
      endif()

      if(BLAS_FOUND)
        message(STATUS "Found BLAS for use with Eigen3")
      endif()
    endif()
  endif()

  set(ADD_EIGEN_SEARCHED_EXTERNAL_LINALG_LIBRARIES TRUE)
endif()

function(target_link_eigen target_name mode)
  # Append the required properties to the passed target
  get_target_property(target_type ${target_name} TYPE)
  if(MKL_FOUND AND SCINE_USE_INTEL_MKL)
    target_link_libraries(${target_name} ${mode} Eigen3::Eigen ${MKL_LIBRARIES} OpenMP::OpenMP_CXX)
    target_include_directories(${target_name} ${mode} ${MKL_INCLUDE_DIRS})
    target_compile_definitions(${target_name} ${mode} EIGEN_USE_MKL_ALL)
  else()
    if(LAPACK_FOUND AND SCINE_USE_LAPACK)
      target_link_libraries(${target_name} ${mode} Eigen3::Eigen ${LAPACK_LIBRARIES})
      target_compile_definitions(${target_name} ${mode} EIGEN_USE_LAPACK EIGEN_USE_BLAS)
    else()
      if(BLAS_FOUND AND SCINE_USE_BLAS)
        target_link_libraries(${target_name} ${mode} Eigen3::Eigen ${BLAS_LIBRARIES})
        target_compile_definitions(${target_name} ${mode} EIGEN_USE_BLAS)
      else()
        target_link_libraries(${target_name} ${mode} Eigen3::Eigen)
      endif()
    endif()
  endif()
endfunction()
