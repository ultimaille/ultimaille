if (BUILD_SHARED_LIBS)
	add_definitions(-DGEO_DYNAMIC_LIBS)
endif ()

add_library(OpenNL_psm OpenNL_psm.c)
target_link_libraries(OpenNL_psm ${CMAKE_DL_LIBS})
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(OpenNL_psm OpenMP::OpenMP_CXX)
endif()
set_target_properties(OpenNL_psm PROPERTIES SOVERSION 1.7)


