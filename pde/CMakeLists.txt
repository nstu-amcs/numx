add_library(pde
  src/sim/sim.c
  src/sim/sim_elm.c
  src/sim/sim_gns.c

  src/msh/msh.c
  src/msh/msh_grd.c
  src/msh/msh_gns.c

  src/fem/fem.c
  src/fem/std/fem.c
  src/fem/std/lin/slv.c
  src/fem/std/lin/asm.c
  src/fem/hmc/fem.c
  src/fem/hmc/lin/slv.c
  src/fem/hmc/lin/asm.c
)

add_library(cgns SHARED IMPORTED)
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/usr/lib/libcgns.so")

target_link_libraries(pde PUBLIC com vec iss cgns)
target_include_directories(pde PUBLIC include/)

add_executable(test_pde 
  test/test.c
  test/ell/test.c
  test/pbc/test.c
  test/hmc/test.c
)
target_link_libraries(test_pde pde munit)
add_test(NAME numx/pde COMMAND test_pde WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

install(DIRECTORY include/numx DESTINATION /usr/share/include)
install(DIRECTORY include/numx DESTINATION /usr/include)

install(FILES cfg/ell.xml DESTINATION /usr/local/share/ElmerGUI/edf)
install(FILES cfg/pbc.xml DESTINATION /usr/local/share/ElmerGUI/edf)
install(FILES cfg/hyp.xml DESTINATION /usr/local/share/ElmerGUI/edf)

add_custom_target(docs-pde-fem-non DEPENDS pde/docs/fem-non.pdf VERBATIM)
add_custom_command(
  OUTPUT pde/docs/fem-non.pdf 
  COMMAND lualatex ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-non.tex 
  DEPENDS ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-non.tex
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_custom_target(docs-pde-fem-pbc DEPENDS pde/docs/fem-pbc.pdf VERBATIM)
add_custom_command(
  OUTPUT pde/docs/fem-pbc.pdf 
  COMMAND lualatex ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-pbc.tex 
  DEPENDS ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-pbc.tex
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_custom_target(docs-pde-fem-hmc DEPENDS pde/docs/fem-hmc.pdf VERBATIM)
add_custom_command(
  OUTPUT pde/docs/fem-hmc.pdf 
  COMMAND lualatex ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-hmc.tex 
  DEPENDS ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-hmc.tex
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
