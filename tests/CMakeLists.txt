cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

Include(FetchContent)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.3.1)

FetchContent_MakeAvailable(Catch2)

FILE(GLOB SRCTEST *.cpp)
add_definitions(-DTEST_INPUT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/input/\")
add_executable(um-test-all ${SRCTEST})
target_link_libraries(um-test-all PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2WithMain)

include(Catch)
catch_discover_tests(um-test-all)

