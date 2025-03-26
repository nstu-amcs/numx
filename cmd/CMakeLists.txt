add_executable(numx
  src/numx.c
  src/pde.c
)

target_link_libraries(numx PRIVATE pde)
target_include_directories(numx PRIVATE include/)

install(TARGETS numx)
