set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)

add_executable(test-permutation test-permutation.cpp)
target_link_libraries(test-permutation m ${CMAKE_DL_LIBS} ultimaille)
