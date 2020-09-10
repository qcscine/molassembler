# Dependencies modeled as conan packages
if(TARGET RingDecomposerLib)
  set(RDL_TARGET "RingDecomposerLib" PARENT_SCOPE)
else()
  add_subdirectory(RingDecomposerLib)
  set(RDL_TARGET "RingDecomposerLibStatic" PARENT_SCOPE)
endif()
