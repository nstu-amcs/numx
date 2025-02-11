file(GLOB SRC src/*.c)

add_library(mtx ${SRC})
target_link_libraries(mtx PUBLIC vec non -lm)
target_include_directories(mtx PUBLIC include/)

add_executable(dmtx_test test/dmtx_test.c)
target_link_libraries(dmtx_test PRIVATE mtx munit)

add_test(NAME dmtx COMMAND dmtx_test)
