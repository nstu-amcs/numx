cmake_minimum_required(VERSION 3.11)
project(numx C)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_C_FLAGS "-std=c23 -Wall -Werror -Wextra -fopenmp")
set(CMAKE_C_STANDARD 23)

# add_compile_definitions(OMP=1)
add_compile_definitions(OMP_NUM_THREADS=1)

add_library(munit munit/munit.c)
target_include_directories(munit PUBLIC munit/)

include(CTest)
enable_testing()

add_subdirectory(com)
add_subdirectory(vec)
add_subdirectory(non)
add_subdirectory(pde)
add_subdirectory(inv)
