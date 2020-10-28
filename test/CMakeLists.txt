OPTION (BIN_INTER "2d interpolation binary" Off)
OPTION (BIN_PGP "2d PGP binary" Off)
OPTION (BIN_PERMU "Permutation test binary " Off)
OPTION (BIN_EMBED "Embeding  binaries" Off)
OPTION (BIN_TEST "Basic tests binary" On)

OPTION (BIN_IO "Basic IO tests" On)





if(BIN_INTER)
	add_subdirectory(2d-interpolation)
endif()
if(BIN_PGP)
	add_subdirectory(pgp2d)
endif()
if(BIN_PERMU)
	add_subdirectory(test-permutation)
endif()
if(BIN_EMBED)
	add_subdirectory(embed)
endif()
if(BIN_TEST)
	add_subdirectory(test)
endif()
if(BIN_IO)
	add_subdirectory(meshIO)
endif()
