cmake_minimum_required(VERSION 2.8)
project(ultimaille-suite)

if(WIN32)
        set(MSVC_CONFIG \$\(Configuration\))
        set(RELATIVE_BIN_DIR bin/${MSVC_CONFIG}/)
        set(RELATIVE_LIB_DIR lib/${MSVC_CONFIG}/)
else()
        set(RELATIVE_BIN_DIR bin/)
        set(RELATIVE_LIB_DIR lib/)
endif()

# Allow the developer to select if Dynamic or Static libraries are built
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" ON)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_LIB_DIR}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/ext)
add_subdirectory(ext)




file(GLOB SOURCES ultimaille/*.cpp)
add_library(ultimaille ${SOURCES})
set_property(TARGET ultimaille PROPERTY CXX_STANDARD 17)
target_link_libraries(ultimaille zlib)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ultimaille OpenMP::OpenMP_CXX)
endif()
set_target_properties(ultimaille PROPERTIES SOVERSION 1.0)



add_subdirectory(test)
