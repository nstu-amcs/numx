add_library(pde
  src/sim.c
  src/msh.c
  src/fem.c

  src/fem_ell.c
  src/fem_pbc.c
  src/fem_hyp.c

  src/msh_grd.c
  src/msh_gns.c
)

target_link_libraries(pde PUBLIC com vec iss)
target_include_directories(pde PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
