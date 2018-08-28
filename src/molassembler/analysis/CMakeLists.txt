project(molassemblerAnalyses)

find_package(Boost REQUIRED COMPONENTS regex program_options system)

file(GLOB analysisFiles "*.cpp")
foreach(filename ${analysisFiles})
  set(targetName "")
  get_filename_component(targetName ${filename} NAME_WE)
  set(targetName "analysis_${targetName}")

  # add the executable and specify compilation and linking options
  add_executable(${targetName} ${filename})

  target_link_libraries(${targetName}
    PRIVATE
      Boost::boost
      Boost::program_options
      Boost::regex
      Boost::system
      molassembler
      gor1
  )
endforeach()
