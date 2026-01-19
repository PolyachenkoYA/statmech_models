//
// Created by ypolyach on 12/1/21.
//

#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

double pown(double x, int n)
{
    double y = x;
    for(int i = 1; i < n; ++i) y *= x;
    return y;
}

int set_f0(double *f, int Nx)
{
    for(int i = 0; i <= Nx; ++i){
        f[i] = sin((((double)i) / Nx) * M_PI);
//        f[i] = exp(-pown(((i * 2.0) / Nx - 1) / 0.2, 2) );
    }

    f[0] = 0;
    f[Nx] = 0;

    return 0;
}

int comp_d2_C(double *f, int Nx, double dx, double* d2)
{
    double dx2 = dx * dx;

    for(int i = 1; i < Nx; ++i){
        d2[i] = (f[i+1] + f[i-1] - f[i] * 2) / dx2;
    }
//    d2[0] = ((f[0] - f[1]) + (f[0] - f[3]) + 4 * (f[2] - f[1])) / dx2;
    d2[0] = (2 * f[0] - 5 * f[1] + 4 * f[2] - f[3]) / dx2;
//    d2[Nx] = ((f[Nx] - f[Nx - 1]) + 4 * (f[Nx - 2] - f[Nx - 1]) + (f[Nx] - f[Nx - 3])) / dx2;
    d2[Nx] = (- f[Nx-3] + 4 * f[Nx-2] - 5 * f[Nx-1] + 2 * f[Nx]) / dx2;

    return 0;
}

double V(double x){ return x * x / 2; }

int get_hf(double* f, double Nx, double dx, double (*V)(double), double* d2, double* hf)
{
    comp_d2_C(f, Nx, dx, d2);
    double L = (dx / 2) * Nx;
    for(int i = 0; i <= Nx; ++i) {
//        hf[i] = ((*V)(i * dx - L) * f[i] - d2[i] / 2 * dx);
        hf[i] = ((*V)(i * dx - L) * f[i] - d2[i] / 2);
//        hf[i] = ((*V)(i * dx - L) * f[i] - d2[i] / 2);
//        hf[i] = ((*V)(i * dx - L) * f[i]);
    }

    return 0;
}

int dot_fncs(double* f1, double* f2, int Nx, double* res)
{
    for(int i = 0; i <= Nx; ++i){
        res[i] = f1[i] * f2[i];
    }
}

double integrate_f12(double* f1, double* f2, int Nx, double dx)   // Simpson
{
//    double I = f1[0] * f2[0] + f1[Nx] * f2[Nx];
//    double I4 = f1[Nx - 1] * f2[Nx - 1];
//    double I2 = 0;
//
//    for(int j = 1; j <= Nx / 2 - 1; ++j){
//        I4 += f1[2 * j - 1] * f2[2 * j - 1];
//        I2 += f1[2 * j] * f2[2 * j];
//    }
//
//    return (I + 2 * I2 + 4 * I4) * (dx / 3);

    double I = (f1[0] * f2[0] + f1[Nx] * f2[Nx]) / 2;
    for(int j = 1; j < Nx; ++j){
        I += f1[j] * f2[j];
    }

    return I * dx;
}

double comp_E(double* f, double Nx, double dx, double (*V)(double), double* d2, double* hf)
{
    get_hf(f, Nx, dx, V, d2, hf);
    double H = integrate_f12(f, hf, Nx, dx);
    double norm = integrate_f12(f, f, Nx, dx);

    return H / norm;
}

int comp_evol_C(double* E, double* f, int Nx, int Nt, double L, double d, int my_seed, int verbose=0)
{
    double dx = 2 * L / Nx;
    double* d2 = (double *) malloc(sizeof (double ) * Nx);
    double* hf = (double *) malloc(sizeof (double ) * Nx);

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, my_seed);

    set_f0(f, Nx);
    E[0] = comp_E(f, Nx, dx, V, d2, hf);

    double f_old;
    int i_move;
    for(int i = 1; i < Nt; ++i){
//        i_move = gsl_rng_uniform_int(rng, Nx + 1);
        i_move = gsl_rng_uniform_int(rng, Nx - 1) + 1;
        f_old = f[i_move];
        f[i_move] += (gsl_rng_uniform(rng) - 0.5) * d;
        E[i] = comp_E(f, Nx, dx, V, d2, hf);
        if(E[i] >= E[i - 1]){
            f[i_move] = f_old;
            E[i] = E[i - 1];
        }
        if(verbose)
            if(i % 1000 == 0)
                printf("progress: %lf\n", (double)(i + 1) / Nt);
    }

    free(d2);
    free(hf);
    gsl_rng_free(rng);
}

pybind11::array_t<double> compute_evolution(int Nx, int Nt, double L, double d, int my_seed, int verbose)
{
    auto Ef = pybind11::array_t<double>(Nt + Nx + 1);   // it will contain both E and M
    pybind11::buffer_info info = Ef.request();
    double *E_ptr = static_cast<double *>(info.ptr);
    double *f_ptr = &(E_ptr[Nt]);

    comp_evol_C(E_ptr, f_ptr, Nx, Nt, L, d, my_seed, verbose);

    return Ef;
}

pybind11::array_t<double> comp_d2(pybind11::array_t<double> f, double dx)
{
    pybind11::buffer_info f_buf = f.request();
    int Nx = f_buf.size - 1;
    pybind11::array_t<double> d2 = pybind11::array_t<double>(Nx + 1);
    pybind11::buffer_info d2_buf = d2.request();
    double *f_ptr = static_cast<double *>(f_buf.ptr);
    double *d2_ptr = static_cast<double *>(d2_buf.ptr);

    comp_d2_C(f_ptr, Nx, dx, d2_ptr);

    return d2;
}
