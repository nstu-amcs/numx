set(SRC
    src/geo.c

    src/sse_fdm.c
)

add_library(pde ${SRC})
target_link_libraries(pde PUBLIC stdx vec iss)
target_include_directories(pde PUBLIC include/)

add_executable(test_geo test/geo/test_geo.c)
target_link_libraries(test_geo PRIVATE pde munit)

add_executable(test_sse test/sse/test_sse.c)
target_link_libraries(test_sse PRIVATE pde munit)

add_test(NAME numx/pde/sse COMMAND test_sse WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/sse)
add_test(NAME numx/pde/geo COMMAND test_geo WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/geo)