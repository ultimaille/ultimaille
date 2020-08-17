set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)

add_executable(embed embed.cpp)
target_link_libraries(embed m ${CMAKE_DL_LIBS} ultimaille)
