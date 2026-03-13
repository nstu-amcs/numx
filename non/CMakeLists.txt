add_library(non 
  src/apx.c
  src/dif.c
  src/opm.c
)

target_link_libraries(non com vec iss)
target_include_directories(non PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
install(DIRECTORY include/numx DESTINATION /usr/include)

add_executable(e_apx_int_cub example/apx/int/cub/main.c)
target_link_libraries(e_apx_int_cub non)

