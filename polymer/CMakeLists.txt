cmake_minimum_required(VERSION 3.20)
project(polymer)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")

set(MATH_LIBS m gsl gslcblas)
set(PY_LIBS crypt pthread dl util rt python3.8)

find_package(PythonLibs REQUIRED)
include_directories(/home/ypolyach/anaconda3/include)
include_directories(/home/ypolyach/anaconda3/include/python3.8)
link_directories(/home/ypolyach/anaconda3/lib)

#add_executable(polymer main.cpp polymer.h polymer.cpp polymer_py.cpp)
add_executable(polymer main.cpp polymer.h polymer.cpp)

target_link_libraries(polymer ${MATH_LIBS} ${PY_LIBS})

find_package(pybind11 REQUIRED)
pybind11_add_module(polymer.so polymer.h polymer.cpp polymer_py.cpp)
target_link_libraries(polymer.so PRIVATE ${MATH_LIBS} ${PY_LIBS})