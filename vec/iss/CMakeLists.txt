add_library(iss
  src/iss.c
  src/diss.c
  src/siss.c
)

target_link_libraries(iss PUBLIC vec mtx dsa)
target_include_directories(iss PUBLIC include/)
