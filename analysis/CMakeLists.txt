# Early return in case analysis binaries are not to be compiled
if(NOT MOLASSEMBLER_ANALYSIS)
  return()
endif()

if(NOT TARGET Boost::regex
  OR NOT TARGET Boost::program_options
  OR NOT TARGET Boost::system
  OR NOT TARGET Boost::filesystem
)
  find_package(Boost REQUIRED
    COMPONENTS regex program_options system filesystem
    QUIET
  )
endif()

if(MOLASSEMBLER_PARALLELIZE AND NOT TARGET OpenMP::OpenMP_CXX)
  find_package(OpenMP QUIET)
endif()

# Build all C++ files in this directory and below
file(GLOB_RECURSE analysisFiles ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(filename ${analysisFiles})
  get_filename_component(targetName ${filename} NAME_WE)
  add_executable(${targetName} ${filename})
  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::program_options
      Boost::regex
      Boost::system
      Boost::filesystem
      ${MOLASSEMBLER_STATIC_TARGET}
      $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )
  target_compile_options(${targetName} PRIVATE ${MOLASSEMBLER_CXX_FLAGS})
endforeach()
