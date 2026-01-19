cmake_minimum_required(VERSION 3.16)
project(test_LANGUAGES_CXX)
set(CMAKE_CXX_STANDARD 14)

set(MATH_LIBS m gsl gslcblas)

add_executable(Izing test.cpp main.cpp)
target_link_libraries(Izing ${MATH_LIBS})

find_package(pybind11 REQUIRED)
pybind11_add_module(test test.cpp py.cpp)
target_link_libraries(test PRIVATE ${MATH_LIBS})
