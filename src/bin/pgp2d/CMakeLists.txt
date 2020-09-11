set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)

add_executable(pgp2d pgp2d.cpp)
target_link_libraries(pgp2d m ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)
