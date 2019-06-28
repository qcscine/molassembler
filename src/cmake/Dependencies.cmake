find_package(Eigen3 3.1.1 REQUIRED)
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

find_package(Boost
  COMPONENTS regex unit_test_framework program_options filesystem system
)
