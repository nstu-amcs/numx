set(SRC
    src/obj.c
)

add_library(geo ${SRC})
target_link_libraries(geo PUBLIC stdx)
target_include_directories(geo PUBLIC include/)

add_executable(test_obj test/test_obj.c)
target_link_libraries(test_obj PRIVATE geo munit)

add_test(NAME test_numx_geo_obj COMMAND test_obj WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test)
