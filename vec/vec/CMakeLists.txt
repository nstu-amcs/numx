add_library(vec 
  src/vec.c
)

target_link_libraries(vec PUBLIC -lm com)
target_include_directories(vec PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)

