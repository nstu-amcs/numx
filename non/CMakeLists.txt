add_library(non 
  src/apx.c
  src/dif.c
  src/non.c
)

target_link_libraries(non com vec dss)
target_include_directories(non PUBLIC include/)
