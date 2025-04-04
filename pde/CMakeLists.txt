add_library(pde
  src/sim.c
  src/sim_elm.c
  src/sim_gns.c

  src/msh.c
  src/msh_grd.c
  src/msh_gns.c

  src/fem.c
  src/fem_ell.c
  src/fem_ell_lin_std.c
  src/fem_ell_lin_non.c
)

add_library(cgns SHARED IMPORTED)
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/usr/lib/libcgns.so")

target_link_libraries(pde PUBLIC com vec iss cgns)
target_include_directories(pde PUBLIC include/)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
