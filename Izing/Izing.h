//
// Created by ypolyach on 10/27/21.
//

#ifndef IZING_IZING_H
#define IZING_IZING_H

#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

int comp_E(int** s, int N, double *E);
int generate_state(int **s, int N, gsl_rng *rng, int mode);
int md(int i, int N);
double get_dE(int **s, int N, int ix, int iy);
int compute_evolution_C(int N, double Temp, int Nt, int my_seed, double *E, double *M);
//void compute_evolution(pybind11::array_t<double> E_py, int N, double Temp, int my_seed=time(nullptr), int verbose=0);
//pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>
pybind11::array_t<double>
        compute_evolution(int N, double Temp, int Nt, int my_seed=time(nullptr), int verbose=0);
int print_E(double *E, int Nt, char prefix=0, char suffix='\n');
int print_S(int **s, int N, char prefix=0);

#endif //IZING_IZING_H
