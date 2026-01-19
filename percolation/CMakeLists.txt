cmake_minimum_required(VERSION 3.19)
project(percolation)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#add_executable(percolation main.cpp percol.hpp percol.cpp percol_v.cpp percol_v.hpp)
#add_executable(percolation main.cpp percol_v.cpp percol_v.hpp)
add_executable(percolation main.cpp percol.cpp percol.hpp)

target_link_libraries(percolation m gsl gslcblas)

#find_package(pybind11 REQUIRED)

#pybind11_add_module(percol prob_pb.cpp prob_for_p.h prob_for_p.cpp)

#target_link_libraries(percol PRIVATE m gsl gslcblas)

