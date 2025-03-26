add_library(com
  src/cmp.c
  src/cut.c
  src/log.c
  src/que.c
)

target_include_directories(com PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
