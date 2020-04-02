if(NOT TARGET Boost::unit_test_framework OR NOT TARGET Boost::regex)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework regex QUIET)
endif()

# Explicitly set all validation files
file(GLOB_RECURSE validationFiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

include(LibType)
determine_lib_type(Boost::unit_test_framework)

# Testing targets
foreach(filename ${validationFiles})
  get_filename_component(targetName ${filename} NAME_WE)
  add_executable(${targetName} ${filename})
  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::regex
      Boost::unit_test_framework
      ${MOLASSEMBLER_STATIC_TARGET}
  )

  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
  if(${LIB_TYPE} STREQUAL "SHARED_LIBRARY")
    target_compile_definitions(${targetName} PRIVATE BOOST_ALL_DYN_LINK)
  endif()
endforeach()
