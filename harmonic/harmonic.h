//
// Created by ypolyach on 12/1/21.
//

#ifndef HARMONIC_HARMONIC_H
#define HARMONIC_HARMONIC_H

#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

double pown(double x, int n);
int set_f0(double *f, int Nx);
int comp_d2_C(double *f, int Nx, double dx, double* d2);
int get_hf(double* f, double Nx, double dx, double (*V)(double), double* d2, double* hf);
int dot_fncs(double* f1, double* f2, int Nx, double* res);
double integrate_f12(double* f1, double* f2, int Nx, double dx);   // Simpson
double comp_E(double* f, double Nx, double dx, double (*V)(double), double* d2, double* hf);
double V(double x);
int mutate_f(double* f, int Nx, double d, gsl_rng* rng);
int comp_evol_C(double* E, double* f, int Nx, int Nt, double L, double d, int my_seed, int verbose=0);

pybind11::array_t<double> compute_evolution(int Nx, int Nt, double L, double d, int my_seed, int verbose);
pybind11::array_t<double> comp_d2(pybind11::array_t<double> f, double dx);

#endif //HARMONIC_HARMONIC_H
