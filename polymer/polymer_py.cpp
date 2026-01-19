//
// Created by ypolyach on 12/1/21.
//

#include <ctime>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "polymer.h"

namespace py = pybind11;

PYBIND11_MODULE(polymer, m)
{
// comp_evol_C(double* E, double* f, int Nx, int Nt, double L, double d, int my_seed, int verbose=0)
// pybind11::array_t<double> compute_evolution(int Nx, int Nt, double L, double d, int my_seed, int verbose)
m.def("compute_evolution", &compute_evolution,
"Fill preallocated array E with 'time'-evolution with Temp for Nt steps; init with my_seed (default=time(NULL))",
py::arg("Nx"),
py::arg("N_timesteps"),
py::arg("my_seed")=time(nullptr),
py::arg("verbose")=0
);
}
