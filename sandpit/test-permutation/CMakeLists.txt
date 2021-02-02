
add_executable(test-permutation test-permutation.cpp)
target_link_libraries(test-permutation  ${CMAKE_DL_LIBS} ultimaille)
IF (NOT WIN32)
  target_link_libraries(test-permutation m)
ENDIF()