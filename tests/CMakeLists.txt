cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-all catch-main.cpp
               test-attributes.cpp
               test-io-geogram.cpp
               test-io-medit.cpp
               test-io-obj.cpp
               test-io-vtk.cpp
               test-numbering.cpp
               test-surface-connectivity.cpp
               test-volume-connectivity.cpp
              )
target_link_libraries(test-all PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

