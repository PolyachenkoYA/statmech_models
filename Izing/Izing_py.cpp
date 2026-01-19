//
// Created by ypolyach on 10/27/21.
//

#include <pybind11/pybind11.h>
#include <ctime>
#include "Izing.h"

namespace py = pybind11;

PYBIND11_MODULE(izing, m)
{
//void compute_evolution(pybind11::array_t<double> E_py, double Temp, int Nt, int my_seed)
    m.def("compute_evolution", &compute_evolution,
        "Fill preallocated array E with 'time'-evolution with Temp for Nt steps; init with my_seed (default=time(NULL))",
        py::arg("grid_size"),
        py::arg("Temp"),
        py::arg("N_steps"),
        py::arg("my_seed")=time(nullptr),
        py::arg("verbose")=0
    );
}

// # double prob_for_p(double p, int box_size, int N_iter, int my_seed=time(NULL), bool verbose=0);
