add_library(msh
  src/umsh/umsh.c
  src/umsh/umsh_grd.c
  src/umsh/umsh_gns.c
  src/umsh/umsh_tel.c
)

target_link_libraries(msh PUBLIC com vec cgns)
target_include_directories(msh PUBLIC include/)

add_executable(umsh
    src/umsh/cmd.c
)

target_link_libraries(umsh PUBLIC msh)
install(TARGETS umsh DESTINATION bin)
