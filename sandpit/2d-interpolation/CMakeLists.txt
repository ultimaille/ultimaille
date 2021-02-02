add_executable(2d-interpolation 2d-interpolation.cpp)
target_link_libraries(2d-interpolation ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)


IF (NOT WIN32)
  target_link_libraries(2d-interpolation m)
ENDIF()