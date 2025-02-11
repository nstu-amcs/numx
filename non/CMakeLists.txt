file(GLOB SRC src/*.c)

add_library(non ${SRC})

target_link_libraries(non vec dss stdx)
target_include_directories(non PUBLIC include/)

add_executable(apx_test test/apx_test.c)

target_link_libraries(apx_test non munit)

add_test(NAME apx COMMAND apx_test)
