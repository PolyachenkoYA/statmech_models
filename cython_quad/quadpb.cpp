#include <pybind11/pybind11.h>
#include "quadratic.hpp"

namespace py = pybind11;

PYBIND11_MODULE(quad, m)
{
    py::class_<Roots>(m, "Roots")
        .def(py::init())
        .def_readwrite("x1", &Roots::x1)
        .def_readwrite("x2", &Roots::x2);

    m.def("quadratic", &quadratic);
}
