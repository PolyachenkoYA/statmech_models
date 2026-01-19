//
// Created by ypolyach on 11/24/21.
//

#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "soft_spheres.h"

int print_state(double *x, double *y, double *z, int N)
{
    for(int i = 0; i < N; ++i){
        printf("%d) (%lf, %lf, %lf)\n", i, x[i], y[i], z[i]);
    }
}

void make_d3(int i, double *xa, double *ya, double *za, double x, double y, double z)
{
    xa[i] = x;
    ya[i] = y;
    za[i] = z;
}

void d3_minus(double* x, double* y, double* z, int i, int j, double* dr)
{
    dr[0] = x[i] - x[j];
    dr[1] = y[i] - y[j];
    dr[2] = z[i] - z[j];
}

double d1_modL(double x, double L)
{
    return (x < 0 ? x + L : (x > L ? x - L : x));
}

void d3_modL(double *dr, double L)
{
    dr[0] = d1_modL(dr[0], L);
    dr[1] = d1_modL(dr[1], L);
    dr[2] = d1_modL(dr[2], L);
}

double d3_dist2(double *dr)
{
    return  dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
}

int powi(int x, int p)
{
    int y = 1;
    for(int i = 0; i < p; ++i) y *= x;
    return y;
}

double powi(double x, int p)
{
    double y = 1.0;
    for(int i = 0; i < p; ++i) y *= x;
    return y;
}

void allocate_r(double *x, double *y, double *z, int N)
{
    x = (double *) malloc(sizeof (double) * N);
    y = (double *) malloc(sizeof (double) * N);
    z = (double *) malloc(sizeof (double) * N);
}

int generate_init(double* x, double* y, double* z, int s, double n, double *L, int *N)
{
    int i, i2, i3;

    double a = pow(4 / n, 1.0 / 3); // L = s * a
//    s = int(pow(N / 4, 1.0 / 3) + 0.5); // it's ok because we check if N is a perfect cube before
    *N = 4 * powi(s, 3);
    *L = a * s;
//    int n_shift = s / 2;
    for(i = 0; i < s; ++i) for(i2 = 0; i2 < s; ++i2) for(i3 = 0; i3 < s; ++i3) { // fcc lattice
                make_d3((i3 + (i2 + i * s) * s) * 4, x, y, z,
                        (i) * a, (i2) * a, (i3) * a);
                make_d3((i3 + (i2 + i * s) * s) * 4 + 1, x, y, z,
                        (i + 0.5) * a, (i2 + 0.5) * a, (i3) * a);
                make_d3((i3 + (i2 + i * s) * s) * 4 + 2, x, y, z,
                        (i + 0.5) * a, (i2) * a, (i3 + 0.5) * a);
                make_d3((i3 + (i2 + i * s) * s) * 4 + 3, x, y, z,
                        (i) * a, (i2 + 0.5) * a, (i3 + 0.5) * a);
            }


    return  0;
}

void compute_total_EP(double* x, double* y, double* z, int N, double L, double Temp, double *U, double *P)
{
    int i, j;
    double r2;
    double _U = 0;
    double _P = 0;
    double dr[3];

    for(i = 0; i < N; ++i){
//        for(j = i + 1; j < N; ++j){
//            d3_minus(x, y, z, i, j, dr);
//            d3_modL(dr, L);
//            r2 = d3_dist2(dr);
//            _U += powi(1.0 / r2, 3);
//        }
        _U += compute_E_ind(x, y, z, N, L, i);
    }

    *U = _U / 2;
//    *U = _U;
    *P = (N * Temp + *U * 6 / 3) / powi(L, 3);
}

double compute_E_ind(double* x, double* y, double* z, int N, double L, int i)
{
    double u = 0;
    double r2, du;
    double dr[3];
    int j;

    for(j = 0; j < N; ++j){
        if(i != j){
            d3_minus(x, y, z, i, j, dr);
            d3_modL(dr, L);
            r2 = d3_dist2(dr);
            du = powi(1.0 / r2, 3);
            u += du;
        }
    }

    return u;
}

void print_Es(double* x, double* y, double* z, int N, double L)
{
    for(int i = 0; i < N; ++i){
        double u = compute_E_ind(x, y, z, N, L, i);
        printf("%d) %lf ", i, u);
    }
    printf("\n");
}

int do_step(double* x, double* y, double* z, int N, double L, double d, double Temp, double *U, double *P, int it, gsl_rng* rng)
{
    int i_move = gsl_rng_uniform_int(rng, N);

//    double u0 = compute_E_ind(x, y, z, N, L, i_move);
    double u0;
    double p0;
    double r_old[3];
    r_old[0] = x[i_move];
    r_old[1] = y[i_move];
    r_old[2] = z[i_move];

    double r_new[3];
    r_new[0] = x[i_move] + (gsl_rng_uniform(rng) - 0.5) * d;
    r_new[1] = y[i_move] + (gsl_rng_uniform(rng) - 0.5) * d;
    r_new[2] = z[i_move] + (gsl_rng_uniform(rng) - 0.5) * d;
    d3_modL(r_new, L);
    x[i_move] = r_new[0];
    y[i_move] = r_new[1];
    z[i_move] = r_new[2];
//    double du = compute_E_ind(x, y, z, N, L, i_move) - u0;
    compute_total_EP(x, y, z, N, L, Temp, &u0, &p0);
    double du = u0 - U[it - 1];

    int accept_step = (du <= 0);
    if(!accept_step) accept_step = (gsl_rng_uniform(rng) < exp(-du / Temp));

//    if(du < -100){
//        printf("it = %d, du = %lf, u0 = %lf\n", it, du, u0);
//        scanf("%s");
//    }
//    if((it == 2182) || (it == 2181) || (it == 2180) || (it == 2692)){
//        printf("it = %d, du = %lf, u0 = %lf, U = %lf\n", it, du, u0, U[it - 1]);
//        print_Es(x, y, z, N, L);
//        print_state(x, y, z, N);
//        getchar();
//    }

    if(accept_step){
//        U[it] = U[it - 1] + du;
        U[it] = u0;
//        P[it] = p0;
        P[it] = (N * Temp + 6 * u0 / 3) / powi(L, 3);
    } else {  // go back
        U[it] = U[it - 1];
        P[it] = P[it - 1];
        x[i_move] = r_old[0];
        y[i_move] = r_old[1];
        z[i_move] = r_old[2];
    }
}

int save_state(char* filename, double* x, double* y, double* z, int N)
{

}

int compute_E_evol_C(int s, double n, int Nt, double d, double Temp, int my_seed, double* U, double* P, int verbose)
{
    int N = 4 * powi(s, 3);
    double L;

    double* x = (double *) malloc(sizeof (double) * N);
    double* y = (double *) malloc(sizeof (double) * N);
    double* z = (double *) malloc(sizeof (double) * N);
    if(verbose) printf("allocated\n");

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, my_seed);

    generate_init(x, y, z, s, n, &L, &N);
    if(verbose) printf("initiated\n");

//    print_state(x, y, z, N);

    compute_total_EP(x, y, z, N, L, Temp, &(U[0]), &(P[0]));
    if(verbose){
        printf("U_0 = %lf; P_0 = %lf\n", U[0], P[0]);
        printf("U_0 done\n");
    }

    int it;
    for(it = 1; it < Nt; ++it){
        do_step(x, y, z, N, L, d, Temp, U, P, it, rng);
        //compute_total_E(x, y, z, N, L, &(U2[it]), &(P[it]));
//        if(verbose){
            if(it % 1000 == 0) {
                printf("progress: %lf\n", double(it + 1) / Nt);
            }
//        }
    }

    free(x);
    free(y);
    free(z);
    gsl_rng_free (rng);

    return 0;
}

pybind11::array_t<double> compute_evolution(int s, double n, int Nt, double d, double Temp, int my_seed, int verbose)
{
//    printf("Nt = %d\n", Nt);

    //    double* E = static_cast<double *>(info.ptr);
    auto EP = pybind11::array_t<double>(Nt * 2);   // it will contain both E and M
    pybind11::buffer_info info = EP.request();
    double *E_ptr = static_cast<double *>(info.ptr);
    double *P_ptr = &(E_ptr[Nt]);

    compute_E_evol_C(s, n, Nt, d, Temp, my_seed, E_ptr, P_ptr, verbose);

    return EP;
}

