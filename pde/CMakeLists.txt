add_library(pde
  src/msh.c
  # src/sim.c
)
target_link_libraries(pde PUBLIC com vec iss)
target_include_directories(pde PUBLIC include/)

