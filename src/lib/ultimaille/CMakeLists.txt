set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_LIB_DIR}/)

file(GLOB SOURCES *.cpp)
add_library(ultimaille SHARED ${SOURCES})
target_link_libraries(ultimaille)
set_target_properties(ultimaille PROPERTIES SOVERSION 1.0)

