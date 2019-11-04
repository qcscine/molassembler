#
# This file is licensed under the 3-clause BSD license.
# Copyright ETH Zurich, Laboratory for Physical Chemistry, Reiher Group.
# See LICENSE.txt for details.
#

# If the target already exists, do nothing
if(NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 REQUIRED)
endif()

# Attempt to find external linalg libraries that accelerate basic calls in
# the following order:
#
# 1. Intel MKL
# 2. LAPACK (brings BLAS too)
# 3. BLAS
#
if(NOT ADD_EIGEN_SEARCHED_EXTERNAL_LINALG_LIBRARIES)
  include(FindMKL)
  if(MKL_FOUND)
    find_package(OpenMP REQUIRED)
    message(STATUS "Found MKL for use with Eigen3")
  else()
    include(FindLAPACK)
    find_package(LAPACK QUIET)
    if(LAPACK_FOUND)
      message(STATUS "Found LAPACK/BLAS for use with Eigen3")
    else()
      include(FindBLAS)
      find_package(BLAS QUIET)
      if(BLAS_FOUND)
        message(STATUS "Found BLAS for use with Eigen3")
      endif()
    endif()
  endif()

  set(ADD_EIGEN_SEARCHED_EXTERNAL_LINALG_LIBRARIES TRUE)
endif()

function(add_eigen target_name)
  # Append the required properties to the passed target
  get_target_property(target_type ${target_name} TYPE)
  if("${target_type}" STREQUAL "INTERFACE_LIBRARY")
    if(MKL_FOUND)
      target_link_libraries(${target_name} INTERFACE Eigen3::Eigen ${MKL_LIBRARIES} OpenMP::OpenMP_CXX)
      target_include_directories(${target_name} INTERFACE ${MKL_INCLUDE_DIRS})
      target_compile_definitions(${target_name} INTERFACE EIGEN_USE_MKL_ALL)
    else()
      if(LAPACK_FOUND)
        target_link_libraries(${target_name} INTERFACE Eigen3::Eigen ${LAPACK_LIBRARIES})
        target_compile_definitions(${target_name} INTERFACE EIGEN_USE_LAPACK EIGEN_USE_BLAS)
      else()
        if(BLAS_FOUND)
          target_link_libraries(${target_name} INTERFACE Eigen3::Eigen ${BLAS_LIBRARIES})
          target_compile_definitions(${target_name} INTERFACE EIGEN_USE_BLAS)
        else()
          target_link_libraries(${target_name} INTERFACE Eigen3::Eigen)
        endif()
      endif()
    endif()
  else()
    if(MKL_FOUND)
      target_link_libraries(${target_name} PUBLIC Eigen3::Eigen ${MKL_LIBRARIES} OpenMP::OpenMP_CXX)
      target_include_directories(${target_name} PRIVATE ${MKL_INCLUDE_DIRS})
      target_compile_definitions(${target_name} PUBLIC EIGEN_USE_MKL_ALL)
    else()
      if(LAPACK_FOUND)
        target_link_libraries(${target_name} PUBLIC Eigen3::Eigen ${LAPACK_LIBRARIES})
        target_compile_definitions(${target_name} PUBLIC EIGEN_USE_LAPACK EIGEN_USE_BLAS)
      else()
        if(BLAS_FOUND)
          target_link_libraries(${target_name} PUBLIC Eigen3::Eigen ${BLAS_LIBRARIES})
          target_compile_definitions(${target_name} PUBLIC EIGEN_USE_BLAS)
        else()
          target_link_libraries(${target_name} PUBLIC Eigen3::Eigen)
        endif()
      endif()
    endif()
  endif()
endfunction()
