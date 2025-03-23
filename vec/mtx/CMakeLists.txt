add_library(mtx 
  src/imtx.c
  src/dmtx.c
  src/smtx.c
  src/jmtx.c
)

target_link_libraries(mtx PUBLIC vec non -lm)
target_include_directories(mtx PUBLIC include/)

