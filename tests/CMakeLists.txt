cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-io test-io.cpp)
target_link_libraries(test-io PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

add_executable(test-tetrahedra test-tetrahedra.cpp)
target_link_libraries(test-tetrahedra PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

