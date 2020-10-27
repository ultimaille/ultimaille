
add_executable(test test.cpp)
target_link_libraries(test m ${CMAKE_DL_LIBS} ultimaille)

add_executable(test-disjoint-set test-disjoint-set.cpp)
target_link_libraries(test-disjoint-set m ${CMAKE_DL_LIBS} ultimaille)

