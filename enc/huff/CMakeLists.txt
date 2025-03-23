add_library(libhuff
  src/huff.c
  src/htree.c
  src/pqueue.c
)

target_include_directories(libhuff PUBLIC src/)

add_executable(huff src/main.c)
target_link_libraries(huff PRIVATE libhuff)

