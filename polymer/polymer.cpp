//
// Created by ypolyach on 12/1/21.
//

#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

int NEIB_INDS[4][2] = {{1,0}, {0,1}, {-1,0}, {0,-1}};

double pown(double x, int n)
{
    double y = x;
    for(int i = 1; i < n; ++i) y *= x;
    return y;
}

void clear_cells(u_int8_t** cells, int Nx)
{
    for(int i = 0; i < Nx; ++i){
        for(int j = 0; j < Nx; ++j){
            cells[i][j] = 0;
        }
    }
}

int pick_next_step(u_int8_t** cells, int Nx, int i, int j, gsl_rng* rng)
{
    int n_free_neib = 0;
    int free_neib_inds[4];
    for(int k = 0; k < 4; ++k){
        if(!cells[i + NEIB_INDS[k][0]][j + NEIB_INDS[k][1]]){
            free_neib_inds[n_free_neib] = k;
            ++n_free_neib;
        }
    }

    if(n_free_neib > 0){
        return free_neib_inds[gsl_rng_uniform_int(rng, n_free_neib)];
    } else {
        return -1;
    }
}

int comp_evol_C(double* R, int Nx, int Nt, int my_seed, int verbose=0)
{
    int L = 2 * Nx + 1;
    u_int8_t** cells = (u_int8_t **) malloc(sizeof (u_int8_t* ) * L);
    for(int i = 0; i < L; ++i){
        cells[i] = (u_int8_t *) malloc(sizeof (u_int8_t ) * L);
    }

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, my_seed);

    int i_t = 0;
    int l, next_step_ind;
//    int n_free_neib = 0;
    int ix, iy;
    int iter = 0;
    while(i_t < Nt){
        clear_cells(cells, L);

        // initial state
        cells[Nx][Nx] = 1;
        cells[Nx][Nx+1] = 1;
        l = 1;
        ix = Nx;
        iy = Nx+1;

        while(l < Nx){
            next_step_ind = pick_next_step(cells, Nx, ix, iy, rng);
            if(next_step_ind > 0){
                ix += NEIB_INDS[next_step_ind][0];
                iy += NEIB_INDS[next_step_ind][1];
                ++l;
                if(verbose){
                    printf("l=%d; (x,y)=(%d,%d)\n", l, ix, iy);
                }
                cells[ix][iy] = 1;
            } else {
                break;
            }
        }

        if(l == Nx){
            R[i_t] = sqrt((ix - Nx) * (ix - Nx) + (iy - Nx) * (iy - Nx));
            ++i_t;
        }

        ++iter;
//        if(verbose)
        if(iter % 1000 == 0)
            printf("progress: %lf\n", (double)(i_t + 1) / Nt);
    }

    free(cells);
    gsl_rng_free(rng);
}

pybind11::array_t<double> compute_evolution(int Nx, int Nt, int my_seed, int verbose)
{
    auto R = pybind11::array_t<double>(Nt);
    pybind11::buffer_info info = R.request();
    double *R_ptr = static_cast<double *>(info.ptr);

    comp_evol_C(R_ptr, Nx, Nt, my_seed, verbose);

    return R;
}
