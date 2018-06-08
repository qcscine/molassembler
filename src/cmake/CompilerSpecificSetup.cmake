set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(UNIX)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} \
      -Wsuggest-attribute=const \
      -Wsuggest-attribute=pure \
      -Wsuggest-final-types \
      -Wsuggest-final-methods"
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-attribute=pure")
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra")
endif()

if(MSVC)
  #Compilation on different cores (faster)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
  #For recursive template depth
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _VARIADIC_MAX=10")
endif()
