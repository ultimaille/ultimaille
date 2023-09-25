add_executable(untangle3d untangle3d.cpp)
target_link_libraries(untangle3d  ${CMAKE_DL_LIBS} ultimaille m)

add_executable(helloworld helloworld.cpp)
target_link_libraries(helloworld  ${CMAKE_DL_LIBS} ultimaille)
IF (NOT WIN32)
	target_link_libraries(helloworld m)
ENDIF()

OPTION (BIN_PGP "2d PGP binary" Off)
OPTION (BIN_EMBED "Embeding  binaries" Off)
OPTION (BIN_TEST "Basic tests binary" Off)

if(BIN_PGP)
	add_subdirectory(pgp2d)
endif()
if(BIN_EMBED)
	add_subdirectory(embed)
endif()
if(BIN_TEST)
	add_subdirectory(test)
endif()

