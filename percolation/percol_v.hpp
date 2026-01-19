//
// Created by ypolyach on 10/15/21.
//

#ifndef PERCOLATION_PERCOL_V_HPP
#define PERCOLATION_PERCOL_V_HPP

#include <ctime>
#include <vector>
#include <gsl/gsl_rng.h>

double prob_for_p(double p, int box_size, int N_iter, int my_seed=time(NULL), bool verbose=0);

int fill_grid_v(std::vector<int>& grid, int N, double p, gsl_rng *rng);
int print_grid_v(std::vector<int>& grid, int N, char prefix=0);
int add_to_cluster_v(std::vector<int>& grid, int N, std::vector<int>& is_checked, std::vector<int>& cluster, int& cluster_size, int pos);
int clear_cluster_v(std::vector<int>& cluster, int& cluster_size);
int is_infinite_cluster_v(std::vector<int>& cluster, int& cluster_size, int N);
int is_infinite_grid_v(std::vector<int>& grid, int N, std::vector<int>& is_checked, std::vector<int>& cluster, int& cluster_size);
int print_cluster_v(std::vector<int>& cluster, int cluster_size, int N, char prefix);

#endif //PERCOLATION_PERCOL_V_HPP
