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

