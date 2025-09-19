#add_library(inv
#  src/inv.c
#)
#
#target_link_libraries(inv PUBLIC vec mtx)
#target_include_directories(inv PUBLIC include/)

add_executable(inv main.c)
target_link_libraries(inv PUBLIC vec mtx dss)
