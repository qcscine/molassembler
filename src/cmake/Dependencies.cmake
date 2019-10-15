find_package(Eigen3 3.1.1 REQUIRED)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

find_package(Boost
  COMPONENTS regex unit_test_framework program_options filesystem system
)

include(FindMKL)
if(MKL_FOUND)
  set(BLAS_FOUND TRUE)
  set(BLAS_LIBRARIES ${MKL_LIBRARIES})
  set(BLAS_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
  set(BLAS_DEFINES EIGEN_USE_MKL EIGEN_USE_BLAS EIGEN_USE_LAPACK)
else()
  include(FindBLAS)
  find_package(BLAS QUIET)
  if(BLAS_FOUND)
    set(BLAS_DEFINES EIGEN_USE_BLAS EIGEN_USE_LAPACK)
  endif()
endif()
