OPTION (BIN_INTER "2d interpolation binary" Off)
OPTION (BIN_PGP "2d PGP binary" On)
OPTION (BIN_PERMU "Permutation test binary " Off)
OPTION (BIN_EMBED "Embeding  binaries" Off)
OPTION (BIN_TEST "Basic tests binary" On)

OPTION (BIN_TEST_IO "Basic IO tests" On)
OPTION (BIN_TEST_VOLUME "Basic volume tests" On)





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
if(BIN_TEST_IO)
	add_subdirectory(test-io)
endif()
if(BIN_TEST_VOLUME)
	add_subdirectory(test-volume)
endif()

