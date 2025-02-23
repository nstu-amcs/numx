cmake_minimum_required(VERSION 3.21)
project(numx C)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_C_FLAGS "-Wall -Werror -Wextra -Wpedantic -std=c23 -fopenmp")
set(CMAKE_C_STANDARD 23)

enable_testing()

add_subdirectory(stdx)
add_subdirectory(vec)
add_subdirectory(non)
add_subdirectory(pde)
# add_subdirectory(opm)
