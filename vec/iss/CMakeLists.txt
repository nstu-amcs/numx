add_library(iss
  src/iss.c
  src/diss.c
  src/siss.c
)

target_link_libraries(iss PUBLIC com vec mtx)
target_include_directories(iss PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
