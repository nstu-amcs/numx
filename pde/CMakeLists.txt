set(SRC
    src/obj.c

    src/sse_fdm.c
)

add_library(pde ${SRC})
target_link_libraries(pde PUBLIC stdx vec iss)
target_include_directories(pde PUBLIC include/)

add_executable(test_obj test/obj/test_obj.c)
target_link_libraries(test_obj PRIVATE pde munit)

add_executable(test_sse test/sse/test_sse.c)
target_link_libraries(test_sse PRIVATE pde munit)

add_test(NAME numx_pde_obj COMMAND test_obj WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/obj)
add_test(NAME numx_pde_sse COMMAND test_sse WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/sse)