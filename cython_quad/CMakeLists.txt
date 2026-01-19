cmake_minimum_required(VERSION 3.16)
project(quadratic_LANGUAGES_CXX)
set(CMAKE_CXX_STANDARD 14)
add_executable(quadratic quadratic.cpp test.cpp quadratic.hpp)
target_link_libraries(quadratic m)

find_package(pybind11 REQUIRED)

pybind11_add_module(quad quadpb.cpp quadratic.cpp quadratic.hpp)

target_link_libraries(quad PRIVATE m)
