cmake_minimum_required(VERSION 2.8)
project(ultimaille-suite)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(RELATIVE_BIN_DIR bin/)
set(RELATIVE_LIB_DIR lib/)

# Allow the developer to select if Dynamic or Static libraries are built
if (WIN32)
    OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
else()
    OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" ON)
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_LIB_DIR}/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${RELATIVE_BIN_DIR}/)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/ext)
add_subdirectory(ext)

if (MSVC)
    # warning level 4 (and all warnings as errors, /WX)
    add_compile_options(/W4)
else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

file(GLOB SOURCES ultimaille/*.cpp ultimaille/*.h ultimaille/algebra/*.cpp ultimaille/algebra/*.h ultimaille/helpers/*.cpp ultimaille/helpers/*.h ultimaille/syntactic-sugar/*.cpp ultimaille/syntactic-sugar/*.h ultimaille/io/*.cpp ultimaille/io/*.h)
add_library(ultimaille ${SOURCES})
#set_property(TARGET ultimaille PROPERTY CXX_STANDARD 17)

target_link_libraries(ultimaille zlib hlbfgs stlbfgs)
if(OpenMP_CXX_FOUND)
    target_link_libraries(ultimaille OpenMP::OpenMP_CXX)
endif()
set_target_properties(ultimaille PROPERTIES SOVERSION 1.0)

add_subdirectory(sandpit)

option(UM_UNIT_TESTS "UM_UNIT_TESTS" OFF)
if (UM_UNIT_TESTS)
    add_subdirectory(tests)
endif()

