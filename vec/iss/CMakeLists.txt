add_library(iss
  src/diss.c
  src/siss.c
)

target_link_libraries(iss PUBLIC com vec mtx)
target_include_directories(iss PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)

# add_executable(test_siss
#   test/siss_test.c
# )
# 
# target_link_libraries(test_siss iss munit)
# add_test(NAME numx/vec/iss/siss COMMAND test_siss WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
