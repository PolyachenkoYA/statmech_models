//
// Created by ypolyach on 11/24/21.
//

#ifndef SOFT_SPHERES_SOFT_SPHERES_H
#define SOFT_SPHERES_SOFT_SPHERES_H

#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

int print_state(double *x, double *y, double *z, int N);
void make_d3(int i, double *xa, double *ya, double *za, double x, double y, double z);
void d3_minus(double* x, double* y, double* z, int i, int j, double* dr);
double d1_modL(double x, double L);
void d3_modL(double *dr, double L);
double d3_dist2(double *dr);
int powi(int x, int p);
double powi(double x, int p);
void allocate_r(double *x, double *y, double *z, int N);
int generate_init(double* x, double* y, double* z, int s, double n, double *L, int *N);
void compute_total_EP(double* x, double* y, double* z, int N, double L, double Temp, double *U, double *P);
double compute_E_ind(double* x, double* y, double* z, int N, double L, int i);
void print_Es(double* x, double* y, double* z, int N, double L);
int do_step(double* x, double* y, double* z, int N, double L, double d, double Temp, double *U, double *P, int it, gsl_rng* rng);
int compute_E_evol_C(int s, double n, int Nt, double d, double Temp, int my_seed, double* U, double* P, int verbose=0);
pybind11::array_t<double> compute_evolution(int s, double n, int Nt, double d, double Temp, int my_seed, int verbose);

#endif //SOFT_SPHERES_SOFT_SPHERES_H
