add_executable(test-io test-io.cpp)
target_link_libraries(test-io ${CMAKE_DL_LIBS} ultimaille)

IF (NOT WIN32)
  target_link_libraries(test-io m)
ENDIF()
