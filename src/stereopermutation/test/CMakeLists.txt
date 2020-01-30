if(NOT SCINE_BUILD_TESTS)
  return()
endif()

if(NOT TARGET Boost::unit_test_framework OR NOT TARGET Boost::program_options)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework program_options)
endif()

add_executable(permutationTests tests.cpp)
target_compile_options(permutationTests PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
if(NOT Boost_USE_STATIC_LIBS)
  target_compile_definitions(permutationTests PRIVATE BOOST_ALL_DYN_LINK)
endif()
target_link_libraries(permutationTests
  PRIVATE
    Boost::boost
    Boost::unit_test_framework
    stereopermutation
)

# add test target
add_test(NAME MASMStereopermutationTests COMMAND permutationTests)

set_target_properties(permutationTests PROPERTIES FOLDER "Tests")
