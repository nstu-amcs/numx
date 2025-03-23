add_library(dss
  src/idss.c
)

target_link_libraries(dss PUBLIC vec mtx)
target_include_directories(dss PUBLIC include/)
