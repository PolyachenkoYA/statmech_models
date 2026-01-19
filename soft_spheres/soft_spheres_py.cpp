//
// Created by ypolyach on 11/24/21.
//
#include <ctime>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "soft_spheres.h"

namespace py = pybind11;

PYBIND11_MODULE(soft_spheres, m)
{
//int compute_E_evol_C(int s, double n, int Nt, double d, double Temp, int my_seed, double* U, double* P, int verbose=0);
//pybind11::array_t<double> compute_evolution(int s, double n, int Nt, double d, double Temp, int my_seed, int verbose);

m.def("compute_evolution", &compute_evolution,
"Fill preallocated array E with 'time'-evolution with Temp for Nt steps; init with my_seed (default=time(NULL))",
py::arg("L/a (number of unit-cells)"),
py::arg("n (concentration)"),
py::arg("N_timesteps"),
py::arg("d (displacement)"),
py::arg("Temp (epsilon units)"),
py::arg("my_seed")=time(nullptr),
py::arg("verbose")=0
);
}
