//
// Created by ypolyach on 10/15/21.
//

#include <pybind11/pybind11.h>
#include "percol.hpp"

namespace py = pybind11;

PYBIND11_MODULE(percol, m)
{
    m.def("prob_for_p", &prob_for_p);
}

