add_library(pde
  src/sim/sim.c
  src/sim/sim_elm.c
  src/sim/sim_cgns.c
  src/sim/sim_tel.c

  src/fem/fem.c
  src/fem/std/std.c
  src/fem/std/lin/lin.c
  src/fem/std/lin/c2d/asm.c
  src/fem/std/lin/c3d/asm.c
  src/fem/hmc/hmc.c
  src/fem/hmc/lin/slv.c
  src/fem/hmc/lin/asm.c
)

target_link_libraries(pde PUBLIC com vec msh iss cgns)
target_include_directories(pde PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
install(DIRECTORY include/numx DESTINATION /usr/include)

add_subdirectory(example)

