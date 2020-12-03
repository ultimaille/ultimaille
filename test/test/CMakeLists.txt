
add_executable(test test.cpp)
target_link_libraries(test  ${CMAKE_DL_LIBS} ultimaille)

add_executable(test-disjoint-set test-disjoint-set.cpp)
target_link_libraries(test-disjoint-set  ${CMAKE_DL_LIBS} ultimaille)

add_executable(test-range test-range.cpp)
set_property(TARGET test-range PROPERTY CXX_STANDARD 17)
target_link_libraries(test-range  ${CMAKE_DL_LIBS} ultimaille)

add_executable(test-lbfgs test-lbfgs.cpp)
target_link_libraries(test-lbfgs  ${CMAKE_DL_LIBS} ultimaille)


IF (NOT WIN32)
	target_link_libraries(test-lbfgs m)
	target_link_libraries(test-range m)
	target_link_libraries(test-disjoint-set m)
	target_link_libraries(test m)
ENDIF()
