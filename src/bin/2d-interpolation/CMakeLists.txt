set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)

add_executable(2d-interpolation 2d-interpolation.cpp)
target_link_libraries(2d-interpolation m ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)


