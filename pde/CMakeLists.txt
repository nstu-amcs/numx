add_library(pde
  src/sim/sim.c
  src/sim/sim_elm.c
  src/sim/sim_gns.c

  src/msh/msh.c
  src/msh/msh_grd.c
  src/msh/msh_gns.c

  src/fem/fem.c
  src/fem/fem_lin_asm.c
  src/fem/fem_lin_slv.c
)

add_library(cgns SHARED IMPORTED)
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/usr/lib/libcgns.so")

target_link_libraries(pde PUBLIC com vec iss cgns)
target_include_directories(pde PUBLIC include/)
target_include_directories(pde PRIVATE src/fem/)

add_executable(test_ell test/ell/test.c)
add_executable(test_non_fpi_def test/non/test_fpi.c)
add_executable(test_non_new_def test/non/test_new.c)
add_executable(test_non_fpi_rlx test/non/test_fpi_rlx.c)
add_executable(test_non_new_rlx test/non/test_new_rlx.c)

target_link_libraries(test_ell pde)
target_link_libraries(test_non_fpi_def pde)
target_link_libraries(test_non_new_def pde)
target_link_libraries(test_non_fpi_rlx pde)
target_link_libraries(test_non_new_rlx pde)

add_test(NAME numx/pde/ell COMMAND test_ell WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME numx/pde/non/fpi/def COMMAND test_non_fpi_def WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME numx/pde/non/new/def COMMAND test_non_new_def WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME numx/pde/non/fpi/rlx COMMAND test_non_fpi_rlx WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME numx/pde/non/new/rlx COMMAND test_non_new_rlx WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

install(DIRECTORY include/numx DESTINATION /usr/share/include)
install(DIRECTORY include/numx DESTINATION /usr/include)

install(FILES mod/ell.xml DESTINATION /usr/local/share/ElmerGUI/edf)
install(FILES mod/pbc.xml DESTINATION /usr/local/share/ElmerGUI/edf)
install(FILES mod/hyp.xml DESTINATION /usr/local/share/ElmerGUI/edf)

add_custom_target(docs-pde-fem-non DEPENDS pde/docs/fem-non.pdf VERBATIM)
add_custom_command(
  OUTPUT pde/docs/fem-non.pdf 
  COMMAND lualatex ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-non.tex 
  DEPENDS ${CMAKE_SOURCE_DIR}/pde/docs/src/fem-non.tex
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
