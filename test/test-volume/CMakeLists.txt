add_executable(test-volume test-volume.cpp)
target_link_libraries(test-volume ${CMAKE_DL_LIBS} ultimaille)

IF (NOT WIN32)
  target_link_libraries(test-volume m)
ENDIF()
