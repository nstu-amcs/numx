add_library(pde
  src/sim.c
  src/sim_elm.c
  src/sim_gns.c

  src/msh.c
  src/msh_grd.c
  src/msh_gns.c

  src/fem.c
  src/fem_lin_asm.c
  src/fem_lin_com.c
  src/fem_lin_ell.c
  src/fem_lin_pbc.c
  src/fem_lin_hyp.c
)

add_library(cgns SHARED IMPORTED)
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/usr/lib/libcgns.so")

target_link_libraries(pde PUBLIC com vec iss cgns)
target_include_directories(pde PUBLIC include/)

add_custom_target(docs-pde-fem-non DEPENDS fem-non.pdf VERBATIM)
add_custom_command(
  OUTPUT fem-non.pdf 
  COMMAND lualatex ${CMAKE_SOURCE_DIR}/pde/docs/fem-non.tex 
  DEPENDS ${CMAKE_SOURCE_DIR}/pde/docs/fem-non.tex
)

install(DIRECTORY include/numx DESTINATION /usr/share/include)

install(FILES model/ell.xml DESTINATION /usr/local/share/ElmerGUI/edf)
install(FILES model/pbc.xml DESTINATION /usr/local/share/ElmerGUI/edf)
install(FILES model/hyp.xml DESTINATION /usr/local/share/ElmerGUI/edf)
