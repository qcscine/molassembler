if(molassembler_FIND_COMPONENTS)
  foreach(component ${molassembler_FIND_COMPONENTS})
    # if molassembler::${component} does not exist, try to import it
    if(NOT TARGET molassembler::${component})
      message(STATUS "Including molassembler component ${component}.")
      include(${CMAKE_CURRENT_LIST_DIR}/${component}/${component}Config.cmake)
    endif()

    # if molassembler::${component} still does not exist, set it to not found
    if(NOT TARGET molassembler::${component})
      set(molassembler_${component}_FOUND 0)
      if(molassembler_FIND_REQUIRED_${component})
        message(FATAL_ERROR "molassembler ${component} not available.")
      endif()
    else()
      set(molassembler_${component}_FOUND 1)
    endif()
  endforeach()
else()
  # By default, include molassembler
  include(${CMAKE_CURRENT_LIST_DIR}/molassembler/molassemblerConfig.cmake)
endif()

