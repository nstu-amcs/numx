set(SRC 
  src/diss.c 
  src/siss.c
)

add_library(iss ${SRC})
target_link_libraries(iss PUBLIC vec mtx)
target_include_directories(iss PUBLIC include/)

add_executable(diss_test test/diss_test.c)
target_link_libraries(diss_test PRIVATE iss munit)

add_test(NAME diss COMMAND diss_test)