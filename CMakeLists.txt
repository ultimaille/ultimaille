cmake_minimum_required(VERSION 3.12...3.26)

project( ultimaille-suite LANGUAGES C CXX )
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release default")
endif()

find_package(OpenMP COMPONENTS CXX)

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

include(FetchContent)
FetchContent_Declare(
  OpenNL_psm
  GIT_REPOSITORY https://github.com/fprotais/geogram.psm.OpenNL
  GIT_TAG        main
)
FetchContent_MakeAvailable(OpenNL_psm)

target_link_libraries(ultimaille zlib stlbfgs OpenNL_psm $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
set_target_properties(ultimaille PROPERTIES SOVERSION 1.0)

add_subdirectory(sandpit)

option(UM_UNIT_TESTS "UM_UNIT_TESTS" OFF)
if (UM_UNIT_TESTS)
    include(CTest)
    enable_testing()
    add_subdirectory(tests)
endif()

file(GENERATE OUTPUT .gitignore CONTENT "*")
