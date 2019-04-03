if(NOT TARGET Boost::unit_test_framework)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework)
endif()
find_package(Eigen3 3.1.1 REQUIRED)

add_executable(templeTests
  ${CMAKE_CURRENT_SOURCE_DIR}/Adaptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cache.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Containers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Functional.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Numeric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Optionals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Variadic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/constexprBTree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/constexprTests.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/various.cpp
)

target_link_libraries(templeTests
  PRIVATE
    Boost::boost
    Boost::unit_test_framework
    temple
)

target_compile_definitions(templeTests PRIVATE BOOST_ALL_DYN_LINK)
target_compile_options(templeTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})

add_test(NAME MASMTempleTests COMMAND templeTests)
set_property(TARGET templeTests PROPERTY FOLDER "Tests")
