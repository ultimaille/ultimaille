
add_executable(pgp2d pgp2d.cpp)
target_link_libraries(pgp2d ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)

add_executable(pgp2d-no-singu pgp2d-no-singu.cpp)
target_link_libraries(pgp2d-no-singu ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)

add_executable(pgp2d-corners pgp2d-corners.cpp)
target_link_libraries(pgp2d-corners ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)

add_executable(cut cut.cpp)
target_link_libraries(cut ${CMAKE_DL_LIBS} ultimaille OpenNL_psm)

IF (NOT WIN32)
  target_link_libraries(pgp2d m)
  target_link_libraries(pgp2d-no-singu m)
  target_link_libraries(pgp2d-corners m)
  target_link_libraries(cut m)
ENDIF()