//
// Created by ypolyach on 12/22/21.
//

#ifndef POLYMER_POLYMER_H
#define POLYMER_POLYMER_H

#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

double pown(double x, int n);
void clear_cells(u_int8_t** cells, int Nx);
int pick_next_step(u_int8_t** cells, int Nx, int i, int j, gsl_rng* rng);
int comp_evol_C(double* R, int Nx, int Nt, int my_seed, int verbose=0);

pybind11::array_t<double> compute_evolution(int Nx, int Nt, int my_seed, int verbose);

#endif //POLYMER_POLYMER_H
