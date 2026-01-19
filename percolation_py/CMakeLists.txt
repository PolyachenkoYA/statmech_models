cmake_minimum_required(VERSION 3.16)
project(percol_LANGUAGES_CXX)
set(CMAKE_CXX_STANDARD 14)
add_executable(percolation percol.cpp percol.hpp test.cpp)
target_link_libraries(percolation m gsl gslcblas)

find_package(pybind11 REQUIRED)

pybind11_add_module(percol percol_pb.cpp percol.cpp percol.hpp)

target_link_libraries(percol PRIVATE m gsl gslcblas)
