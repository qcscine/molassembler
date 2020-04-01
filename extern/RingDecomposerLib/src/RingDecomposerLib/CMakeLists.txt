function(build_library LIBRARY_NAME SPECIFIER)
  add_library("${LIBRARY_NAME}" "${SPECIFIER}"
              RDLutility.c
              RDLhandler.c
              RDLrelation.c
              RDLapsp.c
              RingDecomposerLib.c
              RDLgraph.c
              RDLinfo.c
              RDLcycleFams.c
              RDLstack.c
              RDLtesting.c
              RDLtarjan.c
              RDLdimacs.c
              RDLbitset.c)
  set_property(TARGET "${LIBRARY_NAME}" PROPERTY C_STANDARD 90)
  target_include_directories(
    ${LIBRARY_NAME} PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
  )
  set_target_properties("${LIBRARY_NAME}" PROPERTIES POSITION_INDEPENDENT_CODE True)
  if(NOT WIN32)
    target_link_libraries("${LIBRARY_NAME}" m)
  endif()
  install(
    TARGETS ${LIBRARY_NAME}
    EXPORT RingDecomposerLibTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
  )
endfunction()

build_library("${SHARED_LIBRARY}" "SHARED")
build_library("${STATIC_LIBRARY}" "STATIC")

molassembler_install_component_cmake_files(
  COMPONENT RingDecomposerLib
  EXPORT_NAME RingDecomposerLibTargets
)

install(
  FILES RingDecomposerLib.h
  DESTINATION include
  COMPONENT RingDecomposerLib
)
