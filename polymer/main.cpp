#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "polymer.h"

int main() {
    int Nx = 10;
    int Nt = 100;
    int my_seed = 0;
    int verbose = 1;

    double* R;
    R =  (double*) malloc(sizeof(double) * Nt);
    if(verbose) printf("allocated\n");

    comp_evol_C(R, Nx, Nt, my_seed);
    if(verbose) printf("done\n");

    free(R);
    if(verbose) printf("freed\n");

    return 0;
}
