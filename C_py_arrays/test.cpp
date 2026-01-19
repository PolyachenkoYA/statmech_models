#include <cmath>
#include <ctime>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

int comp_E(int** s, int N, double *E);
int generate_state(int **s, int N, gsl_rng *rng, int mode);
int md(int i, int N);
double get_dE(int **s, int N, int ix, int iy);
int compute_evolution_C(int N, double Temp, int Nt, int my_seed, double *E);
//void compute_evolution(py::array_t<double> E_py, int N, double Temp, int my_seed=time(nullptr), int verbose=0);
//py::array_t<double, py::array::c_style | py::array::forcecast>
py::array_t<double>
        compute_evolution(int N, double Temp, int Nt, int my_seed=time(nullptr), int verbose=0);
int print_E(double *E, int Nt, char prefix=0, char suffix='\n');

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
    _E += s[N-1][N-1] * (s[0][N-1] + s[N-1][0]);

    *E = -_E;

    return 0;
}

int generate_state(int **s, int N, gsl_rng *rng, int mode = 1)
{
    int i, j;

    if(mode == 0){ // random
        for(i = 0; i < N-1; ++i){
            for(j = 0; j < N-1; ++j) {
                s[i][j] = gsl_rng_uniform_int(rng, 2) * 2 - 1;
            }
        }
    } else if(mode == 1) { // ordered
        for(i = 0; i < N-1; ++i){
            for(j = 0; j < N-1; ++j) {
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

int compute_evolution_C(int N, double Temp, int Nt, int my_seed, double *E)
{
    int i;
    //double *E = (double*) malloc(sizeof(double) * Nt);

    print_E(E, Nt, 'c');

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

    int ix, iy;
    double dE;
    int do_flip;
    for(i = 1; i < Nt; ++i){
        ix = gsl_rng_uniform_int(rng, N);
        iy = gsl_rng_uniform_int(rng, N);

        dE = get_dE(s, N, ix, iy);
        do_flip = (dE <= 0 ? 1 : (gsl_rng_uniform(rng) < exp(- dE / Temp) ? 1 : 0));
        if(do_flip){
            s[ix][iy] *= -1;
            E[i] = E[i-1] + dE;
        } else {
            E[i] = E[i-1];
        }
    }

    for(i = 0; i < N; ++i) free(s[i]);
    free(s);
    gsl_rng_free (rng);

    print_E(E, Nt, 'C');

    return 0;
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


