cmake_minimum_required(VERSION 3.11)
project(numx C)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_C_FLAGS "-std=gnu2x -Wall -Wextra -fopenmp  -O2")
set(CMAKE_C_STANDARD 23)

# add_library(munit munit/munit.c)
# target_include_directories(munit PUBLIC munit/)

add_library(cgns SHARED IMPORTED)
set_property(TARGET cgns PROPERTY IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/libcgns.so")
include_directories(SYSTEM "/usr/include")

include(CTest)
enable_testing()

add_subdirectory(com)
add_subdirectory(vec)
add_subdirectory(non)
add_subdirectory(msh)
add_subdirectory(pde)
add_subdirectory(inv)
