cmake_minimum_required(VERSION 2.8)
project(ultimaille-test)

set(ULTIMAILLE_TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR})

# Allow the developer to select if Dynamic or Static libraries are built
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" ON)

add_subdirectory(src)
