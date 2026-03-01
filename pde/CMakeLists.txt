add_library(pde
  src/sim/sim.c
  src/sim/sim_elm.c
  src/sim/sim_gns.c
  src/sim/sim_tel.c

  src/fem/fem.c
  src/fem/std/fem.c
  src/fem/std/lin/3d/slv.c
  src/fem/std/lin/3d/asm.c
  src/fem/hmc/fem.c
  src/fem/hmc/lin/slv.c
  src/fem/hmc/lin/asm.c
)

target_link_libraries(pde PUBLIC com vec msh iss cgns)
target_include_directories(pde PUBLIC include/)

add_executable(e_emf example/emf/main.c)
target_link_libraries(e_emf pde)
add_test(NAME numx/example/emf COMMAND e_emf)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
install(DIRECTORY include/numx DESTINATION /usr/include)
