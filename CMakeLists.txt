set(Boost_USE_STATIC_LIBS ON)

find_package (Boost COMPONENTS unit_test_framework REQUIRED)
find_package (Eigen3 REQUIRED)

include_directories(SYSTEM
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

if( UNIX )
  add_definitions(-fPIC)
endif()

add_library(molecule STATIC ${SOURCE_all} ${HEADERS_all})

set_property(TARGET molecule PROPERTY CXX_STANDARD 14)

target_link_libraries(molecule ${Boost_UNIT_TEST_FRAMEWORK})
