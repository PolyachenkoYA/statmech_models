//
// Created by ypolyach on 10/27/21.
//

#include <gsl/gsl_rng.h>
#include <cmath>
#include <ctime>

#include "Izing.h"

int comp_M(int **s, int N, double *M)
{
    int i, j;
    double _M = 0;
    for(i = 0; i < N; ++i){
        for(j = 0; j < N; ++j){
            _M += s[i][j];
        }
    }

    *M = _M;
}

int comp_E(int** s, int N, double *E)
{
    int i, j;
    double _E = 0;
    for(i = 0; i < N-1; ++i){
        for(j = 0; j < N-1; ++j){
            _E += s[i][j] * (s[i+1][j] + s[i][j+1]);
        }
        _E += s[i][N-1] * (s[i+1][N-1] + s[i][0]);
    }
    for(j = 0; j < N-1; ++j){
        _E += s[N-1][j] * (s[0][j] + s[N-1][j+1]);
    }
    _E += s[N-1][N-1] * (s[0][N-1] + s[N-1][0]);

    *E = -_E;

    return 0;
}

int generate_state(int **s, int N, gsl_rng *rng, int mode = 1)
{
    int i, j;

    if(mode == 0){ // random
        for(i = 0; i < N; ++i){
            for(j = 0; j < N; ++j) {
                s[i][j] = gsl_rng_uniform_int(rng, 2) * 2 - 1;
            }
        }
    } else if(mode == 1) { // ordered
        for(i = 0; i < N; ++i){
            for(j = 0; j < N; ++j) {
                s[i][j] = 1;
            }
        }
    }
}

int md(int i, int N)
{
    return i >= 0 ? (i < N ? i : 0) : (N-1);   // i mod N
}

double get_dE(int **s, int N, int ix, int iy)
{
    return 2 * s[ix][iy] * (s[md(ix + 1, N)][iy]
                            + s[ix][md(iy + 1, N)]
                            + s[md(ix - 1, N)][iy]
                            + s[ix][md(iy - 1, N)]);
}

int compute_evolution_C(int N, double Temp, int Nt, int my_seed, double *E, double *M)
{
    int i;
    //double *E = (double*) malloc(sizeof(double) * Nt);

//    E[0] = 123;
//    print_E(E, Nt, 'c');
//    printf("N = %d\n", N);

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, my_seed);

    int **s = (int**) malloc(sizeof(int*) * N);
    for(i = 0; i < N; ++i){
        s[i] = (int*) malloc(sizeof(int) * N);
    }

    generate_state(s, N, rng, 1);
    comp_E(s, N, &(E[0])); // &(E[0]) == E
    comp_M(s, N, &(M[0]));

//    print_E(E, Nt, '0');

    int ix, iy;
    double dE;
    int do_flip;
    for(i = 1; i < Nt; ++i){
        ix = gsl_rng_uniform_int(rng, N);
        iy = gsl_rng_uniform_int(rng, N);

        dE = get_dE(s, N, ix, iy);
        do_flip = (dE <= 0 ? 1 : (gsl_rng_uniform(rng) < exp(- dE / Temp) ? 1 : 0));
//        printf("%d, %lf\n", do_flip, dE);
//        print_S(s, N, '0' + (i % 10));
        if(do_flip){
            M[i] = M[i-1] - 2 * s[ix][iy];
            E[i] = E[i-1] + dE;
            s[ix][iy] *= -1;
        } else {
            E[i] = E[i-1];
            M[i] = M[i-1];
        }

//        print_E(E, Nt, '0' + (i % 10));
    }

    for(i = 0; i < N; ++i) free(s[i]);
    free(s);
    gsl_rng_free (rng);

//    print_E(E, Nt, 'C');

    return 0;
}

pybind11::array_t<double> compute_evolution(int N, double Temp, int Nt, int my_seed, int verbose)
{
//    printf("Nt = %d\n", Nt);

    //    double* E = static_cast<double *>(info.ptr);
    auto EM = pybind11::array_t<double>(Nt * 2);   // it will contain both E and M
    pybind11::buffer_info info = EM.request();
    double *E_ptr = static_cast<double *>(info.ptr);
    double *M_ptr = &(E_ptr[Nt]);

//    print_E(E_ptr, Nt, 'p');

    compute_evolution_C(N, Temp, Nt, my_seed, E_ptr, M_ptr);

//    print_E(E_ptr, Nt, 'P');

    return EM;
}

int print_E(double *E, int Nt, char prefix, char suffix)
{
    if(prefix > 0){
        printf("%c\n", prefix);
    }

    for(int i = 0; i < Nt; ++i){
        printf("%lf ", E[i]);
    }

    if(suffix > 0){
        printf("%c", suffix);
    }
}

int print_S(int **s, int N, char prefix)
{
    int i, j;

    if(prefix > 0){
        printf("%c\n", prefix);
    }

    for(i = 0; i < N; ++i){
        for(j = 0; j < N; ++j){
            printf("%2d", s[i][j]);
        }
        printf("\n");
    }
}
