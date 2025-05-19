add_library(non 
  src/apx.c
  src/dif.c
  src/opm.c
)

target_link_libraries(non com vec dss)
target_include_directories(non PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
install(DIRECTORY include/numx DESTINATION /usr/include)
