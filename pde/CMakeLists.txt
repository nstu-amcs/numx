set(SRC
    src/geo.c
    src/cnd.c

    src/sse_ctx.c
    src/sse_fdm.c
    src/sse_fem.c
)

add_library(pde ${SRC})
target_link_libraries(pde PUBLIC stdx vec iss)
target_include_directories(pde PUBLIC include/)

add_executable(test_geo test/geo/test_geo.c)
target_link_libraries(test_geo PRIVATE pde munit)

add_executable(test_sse_fdm test/sse/test_sse_fdm.c)
target_link_libraries(test_sse_fdm PRIVATE pde munit)

add_test(NAME numx/pde/geo COMMAND test_geo WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/geo)
add_test(NAME numx/pde/sse/fdm COMMAND test_sse_fdm WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/sse)