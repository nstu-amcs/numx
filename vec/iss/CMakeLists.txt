set(SRC 
  src/iss.c

  src/diss.c 
  src/siss.c
)

add_library(iss ${SRC})
target_link_libraries(iss PUBLIC vec mtx stdx)
target_include_directories(iss PUBLIC include/)

add_executable(test_diss test/test_diss.c)
target_link_libraries(test_diss PRIVATE iss munit)

add_test(NAME numx/vec/iss/diss COMMAND test_diss)