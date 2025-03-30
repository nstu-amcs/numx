add_library(pde
  src/sim.c
  src/msh.c
  src/fem.c

  src/sim_elm.c

  src/msh_grd.c
  src/msh_gns.c

  src/fem_ell.c
  src/fem_pbc.c
  src/fem_hyp.c
)

add_library(cgns SHARED IMPORTED)
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/usr/lib/libcgns.so")

target_link_libraries(pde PUBLIC com vec iss cgns)
target_include_directories(pde PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
