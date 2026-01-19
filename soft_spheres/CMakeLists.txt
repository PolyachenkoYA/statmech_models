cmake_minimum_required(VERSION 3.20)
project(soft_spheres)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(MATH_LIBS m gsl gslcblas)
set(PY_LIBS crypt pthread dl util rt python3.8)

find_package(PythonLibs REQUIRED)
include_directories(/home/ypolyach/anaconda3/include)
include_directories(/home/ypolyach/anaconda3/include/python3.8)
link_directories(/home/ypolyach/anaconda3/lib)

add_executable(soft_spheres main.cpp soft_spheres.cpp soft_spheres.h soft_spheres_py.cpp)
target_link_libraries(soft_spheres ${MATH_LIBS} ${PY_LIBS})

find_package(pybind11 REQUIRED)
pybind11_add_module(soft_spheres.so soft_spheres.h soft_spheres.cpp soft_spheres_py.cpp)
target_link_libraries(soft_spheres.so PRIVATE ${MATH_LIBS} ${PY_LIBS})


