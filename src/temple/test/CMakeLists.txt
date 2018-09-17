project(templeTestProject)

find_package(Boost REQUIRED COMPONENTS unit_test_framework)
find_package(Eigen3 3.1.1 REQUIRED)

add_executable(templeTests
  ${CMAKE_CURRENT_SOURCE_DIR}/Adaptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Cache.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Containers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Functional.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Numeric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Optionals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/constexprBTree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/constexprTests.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/various.cpp
)

if(
  "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"
  OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"
)
  target_compile_options(templeTests
    PRIVATE -fconstexpr-steps=1999999999 -fconstexpr-backtrace-limit=0
  )
endif()

target_link_libraries(templeTests
  PRIVATE
    Boost::boost
    Boost::unit_test_framework
    temple
)

add_test(NAME templeTests COMMAND templeTests)
set_property(TARGET templeTests PROPERTY FOLDER "Tests")
