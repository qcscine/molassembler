if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wpedantic -Wextra")
endif()

if(MSVC)
  #Compilation on different cores (faster)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
  #For recursive template depth
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _VARIADIC_MAX=10")
endif()
