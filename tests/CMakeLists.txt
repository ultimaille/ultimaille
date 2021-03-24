cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v2.13.1)

FetchContent_MakeAvailable(Catch2)

add_executable(test-io-geogram test-io-geogram.cpp)
target_link_libraries(test-io-geogram PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

add_executable(test-io-obj test-io-obj.cpp)
target_link_libraries(test-io-obj PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

add_executable(test-io-medit test-io-medit.cpp)
target_link_libraries(test-io-medit PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

add_executable(test-io-vtk test-io-vtk.cpp)
target_link_libraries(test-io-vtk PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

add_executable(test-numbering test-numbering.cpp)
target_link_libraries(test-numbering PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

add_executable(test-volume-connectivity test-volume-connectivity.cpp)
target_link_libraries(test-volume-connectivity PRIVATE ${CMAKE_DL_LIBS} ultimaille Catch2::Catch2)

