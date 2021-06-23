cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
include(CTest)
enable_testing()

Include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

file(GLOB TEST_SOURCES *.cpp)
add_definitions(-DTEST_INPUT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/input/\")
add_executable(um-test-all ${TEST_SOURCES})
target_link_libraries(um-test-all PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

FetchContent_GetProperties(Catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/contrib)
include(Catch)

catch_discover_tests(um-test-all)
