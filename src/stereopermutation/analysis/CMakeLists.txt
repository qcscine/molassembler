project(analysis)

if(NOT TARGET Boost::program_options)
  find_package(Boost REQUIRED COMPONENTS program_options)
endif()

add_executable(stereopermutationProperties properties.cpp)
target_compile_options(stereopermutationProperties PUBLIC ${MOLASSEMBLER_CXX_FLAGS})
target_link_libraries(stereopermutationProperties
  PRIVATE Boost::program_options stereopermutation
)

set_property(TARGET stereopermutationProperties PROPERTY FOLDER "Analysis")
