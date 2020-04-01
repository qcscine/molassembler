# Early return in case validation tests are not to be compiled
if(NOT MOLASSEMBLER_VALIDATION_TESTS)
  return()
endif()

if(NOT TARGET Boost::unit_test_framework OR NOT TARGET Boost::regex)
  find_package(Boost REQUIRED COMPONENTS unit_test_framework regex QUIET)
endif()

# Explicitly set all test files
set(testFiles
  ${CMAKE_CURRENT_SOURCE_DIR}/GraphAlgorithms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RankingValidation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/StrainedMolecules.cpp
)

# Testing targets
foreach(filename ${testFiles})
  # get a target name from the filename
  set(baseName "")
  get_filename_component(baseName ${filename} NAME_WE)

  set(targetName "validation_${baseName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::regex
      Boost::unit_test_framework
      gor1
      ${MOLASSEMBLER_STATIC_TARGET}
      stereopermutation
  )

  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
endforeach()
