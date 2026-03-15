add_executable(e_pde_mf_non main.c)
target_link_libraries(e_pde_mf_non pde)
add_test(NAME numx/example/pde/mf/non COMMAND e_pde_mf_non)
