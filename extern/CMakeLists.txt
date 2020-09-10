# Dependencies modeled as conan packages
if(TARGET RingDecomposerLib)
  set(RDL_TARGET "RingDecomposerLib" PARENT_SCOPE)
else()
  add_subdirectory(RingDecomposerLib)
  set(RDL_TARGET "RingDecomposerLibStatic" PARENT_SCOPE)
endif()

if(NOT TARGET nauty)
  add_subdirectory(nauty)
endif()

# Single-header dependencies
add_subdirectory(outcome)
add_subdirectory(json)

add_subdirectory(licenses)
