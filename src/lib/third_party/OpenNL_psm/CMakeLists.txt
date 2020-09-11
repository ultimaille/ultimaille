if (BUILD_SHARED_LIBS)
	add_definitions(-DGEO_DYNAMIC_LIBS)
endif ()

add_library(OpenNL_psm OpenNL_psm.c)
set_target_properties(OpenNL_psm PROPERTIES SOVERSION 1.7)


