#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "harmonic.h"

int main() {
    int Nx = 100 * 2;
    int Nt = 100;
    double L = 1.0;
    int my_seed = 0;
    double d = 0.1;
    int verbose = 1;

    double* f;
    double* E;
    f = (double*) malloc(sizeof(double) * (Nx + 1));
    E =  (double*) malloc(sizeof(double) * Nt);
    if(verbose) printf("allocated\n");

    comp_evol_C(E, f, Nx, Nt, L, d, my_seed);
    if(verbose) printf("done\n");

    free(f);
    free(E);
    if(verbose) printf("freed\n");

    return 0;
}
