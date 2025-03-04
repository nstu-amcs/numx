set(SRC src/vec.c)

add_library(vec ${SRC} ${HDR})

target_link_libraries(vec PUBLIC -lm)
target_include_directories(vec PUBLIC include/)

add_executable(vec_test test/vec_test.c)

target_link_libraries(vec_test vec munit)

add_test(NAME numx/vec/vec COMMAND vec_test)
