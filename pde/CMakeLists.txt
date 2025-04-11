add_library(pde
  src/sim.c
  src/sim_elm.c
  src/sim_gns.c

  src/msh.c
  src/msh_grd.c
  src/msh_gns.c

  src/fem.c
  src/fem_lin_com.c
  src/fem_lin_ell.c
  src/fem_lin_pbc.c
)

add_library(cgns SHARED IMPORTED)
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/usr/lib/libcgns.so")

target_link_libraries(pde PUBLIC com vec iss cgns)
target_include_directories(pde PUBLIC include/)

# add_executable(case4 test/case4/test.c)
# target_link_libraries(case4 PRIVATE pde)
# add_test(NAME numx/pde/case4 COMMAND case4 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/pde/test/case4)

install(DIRECTORY include/numx DESTINATION /usr/share/include)
