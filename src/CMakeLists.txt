# RELATIVE_BIN_DIR
# RELATIVE_LIB_DIR
#  Directories where libraries and executables are generated, relative to
# source directory. These variables are to be used as an input (e.g., link
# path, or path to find GOMGEN executable)
#

if(WIN32)
        set(MSVC_CONFIG \$\(Configuration\))
        set(RELATIVE_BIN_DIR bin/${MSVC_CONFIG}/)
        set(RELATIVE_LIB_DIR lib/${MSVC_CONFIG}/)
else()
        set(RELATIVE_BIN_DIR bin/)
        set(RELATIVE_LIB_DIR lib/)
endif()

#
# defines compilation flags
#

include(${ULTIMAILLE_TEST_SOURCE_DIR}/cmake/cxx_compiler_flag.cmake)

enable_cxx_compiler_flag_if_supported("-Wall")
enable_cxx_compiler_flag_if_supported("-Wextra")
enable_cxx_compiler_flag_if_supported("-pedantic")
enable_cxx_compiler_flag_if_supported("-std=c++17")
enable_cxx_compiler_flag_if_supported("-O3")
enable_cxx_compiler_flag_if_supported("-g")
enable_cxx_compiler_flag_if_supported("-fopenmp")
enable_cxx_compiler_flag_if_supported("-mieee-fp")
enable_cxx_compiler_flag_if_supported("-ffp-contract=off")

#
#
#

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib ${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party)

add_subdirectory(lib)
add_subdirectory(bin)

