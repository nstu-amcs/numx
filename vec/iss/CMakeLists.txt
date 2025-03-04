set(SRC 
  src/diss.c 
  src/siss.c
)

add_library(iss ${SRC})
target_link_libraries(iss PUBLIC vec mtx stdx)
target_include_directories(iss PUBLIC include/)

add_executable(diss_test test/diss_test.c)
target_link_libraries(diss_test PRIVATE iss munit)

add_test(NAME numx/vec/iss/diss COMMAND diss_test)