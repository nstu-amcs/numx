add_executable(e_pde_mf_lin main.c)
target_link_libraries(e_pde_mf_lin pde)
add_test(NAME numx/example/pde/mf COMMAND e_pde_mf_lin)
