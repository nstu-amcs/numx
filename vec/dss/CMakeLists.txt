set(SRC 
  src/idss.c
)

add_library(dss ${SRC})

target_link_libraries(dss PUBLIC vec mtx)
target_include_directories(dss PUBLIC include/)

add_executable(test_imp test/test_imp.c)

target_link_libraries(test_imp dss munit)

add_test(NAME numx/vec/dss/idss COMMAND test_imp)
