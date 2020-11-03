add_executable(pgp2d pgp2d.cpp)
target_link_libraries(pgp2d ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)

IF (NOT WIN32)
  target_link_libraries(pgp2d m)
ENDIF()
