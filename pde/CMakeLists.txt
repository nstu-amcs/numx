add_library(pde
  src/grd.c
  src/sim.c

  src/fem.c
  src/fem_ell.c
  src/fem_pbc.c
)

target_link_libraries(pde PUBLIC com vec iss)
target_include_directories(pde PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
