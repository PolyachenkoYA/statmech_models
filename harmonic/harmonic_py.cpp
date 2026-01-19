//
// Created by ypolyach on 12/1/21.
//

//
// Created by ypolyach on 11/24/21.
//
#include <ctime>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "harmonic.h"

namespace py = pybind11;

PYBIND11_MODULE(harmonic, m)
{
// comp_evol_C(double* E, double* f, int Nx, int Nt, double L, double d, int my_seed, int verbose=0)
// pybind11::array_t<double> compute_evolution(int Nx, int Nt, double L, double d, int my_seed, int verbose)
m.def("compute_evolution", &compute_evolution,
"Fill preallocated array E with 'time'-evolution with Temp for Nt steps; init with my_seed (default=time(NULL))",
py::arg("Nx (Nx = L/dx * 2"),
py::arg("N_timesteps"),
py::arg("L (x in [-L...L]))"),
py::arg("d (df in [-d/2; d/2])"),
py::arg("my_seed")=time(nullptr),
py::arg("verbose")=0
);

// pybind11::array_t<double> comp_d2(pybind11::array_t<double> f, double dx)
m.def("comp_d2", &comp_d2,
py::arg("f"),
py::arg("dx")
);
}
