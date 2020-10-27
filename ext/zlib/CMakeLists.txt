if (BUILD_SHARED_LIBS)
	add_definitions(-DZLIB_DLL)
endif ()

file(GLOB SOURCES *.c)
add_library(zlib ${SOURCES})
set_target_properties(zlib PROPERTIES SOVERSION 1.2)
