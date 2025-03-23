add_library(com
  src/cut.c
  src/log.c
  src/que.c
)

target_include_directories(com PUBLIC include/)
