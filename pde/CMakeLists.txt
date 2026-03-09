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

add_executable(e_p1 example/p1/main.c)
target_link_libraries(e_p1 pde)
add_test(NAME numx/example/p1 COMMAND e_p1)

add_executable(e_p2_u example/p2-u/main.c)
target_link_libraries(e_p2_u pde)
add_test(NAME numx/example/p2-u COMMAND e_p2_u)

add_executable(e_p2_n example/p2-n/main.c)
target_link_libraries(e_p2_n pde)
add_test(NAME numx/example/p2-n COMMAND e_p2_n)

add_executable(e_donut example/donut/main.c)
target_link_libraries(e_donut pde)
add_test(NAME numx/example/donut COMMAND e_donut)

add_executable(e_emf example/emf/main.c)
target_link_libraries(e_emf pde)
add_test(NAME numx/example/emf COMMAND e_emf)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
install(DIRECTORY include/numx DESTINATION /usr/include)
