#include <pybind11/pybind11.h>
#include <time.h>
#include "percol.hpp"

namespace py = pybind11;

PYBIND11_MODULE(percol, m)
{
    m.def("prob_for_p", &prob_for_p,
		  "Percolation probability esimation for given cell-fill probability, box size and N_iter",
		   py::arg("p"),
		   py::arg("box_size"),
		   py::arg("N_iter"),
		   py::arg("my_seed")=time(NULL),
		   py::arg("verbose")=0
 		);
}

// # double prob_for_p(double p, int box_size, int N_iter, int my_seed=time(NULL), bool verbose=0);
