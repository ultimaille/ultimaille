add_executable(test-io-volume test-io-volume.cpp)
target_link_libraries(test-io-volume ${CMAKE_DL_LIBS} ultimaille)

IF (NOT WIN32)
  target_link_libraries(test-io-volume m)
ENDIF()