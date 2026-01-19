//
// Created by ypolyach on 10/15/21.
//

#ifndef PERCOLATION_PERCOL_HPP
#define PERCOLATION_PERCOL_HPP

#include <time.h>
#include <gsl/gsl_rng.h>

double prob_for_p(double p, int box_size, int N_iter, int my_seed=time(NULL), bool verbose=0);

int fill_grid(int* grid, int N, double p, gsl_rng *rng);
int print_grid(int* grid, int N, char prefix=0);
int add_to_cluster(int* grid, int N, int* is_checked, int* cluster, int* cluster_size, int pos);
int clear_cluster(int* cluster, int *cluster_size);
int is_infinite_cluster(int* cluster, int* cluster_size, int N);
int is_infinite_grid(int* grid, int N, int* is_checked, int* cluster, int* cluster_size);
int print_cluster(int* cluster, int cluster_size, int N, char prefix);
int copy_arr(int* from, int* to, int N);

#endif //PERCOLATION_PERCOL_HPP
