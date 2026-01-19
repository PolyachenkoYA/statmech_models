//
// Created by ypolyach on 10/15/21.
//

#include <cstdio>
#include <gsl/gsl_rng.h>
#include <ctime>
#include <vector>

#include "percol_v.hpp"

int fill_grid_v2(std::vector<int>& grid, int N, double p)
{
    int N2 = N * N;
    int i;

    for(i = 0; i < N2; ++i){
//        grid[i] = (gsl_rng_uniform(rng) < p);
        grid[i] = ((double)rand() / RAND_MAX < p);
    }

    return 0;
}

double prob_for_p(double p, int box_size, int N_iter, int my_seed, bool verbose)
{
    int box_size2 = box_size * box_size;

    std::vector<int> grid_v;
    grid_v.resize(box_size2);
    std::vector<int> is_checked_v;
    is_checked_v.resize(box_size2);
    std::vector<int> cluster_v;
    cluster_v.resize(box_size2);

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, my_seed);
//    srand(my_seed);

    int i;
    int cluster_size = box_size * box_size;   // initial N2 so it's all initialized to -1
    int res = 0;
    int is_inf_gr;

    for(i = 0; i < N_iter; ++i){
        fill_grid_v(grid_v, box_size, p, rng);
//        fill_grid_v2(grid_v, box_size, p);

        if(verbose) {
            print_grid_v(grid_v, box_size);
        }

        is_inf_gr = is_infinite_grid_v(grid_v, box_size, is_checked_v, cluster_v, cluster_size);
        res += is_inf_gr;

        if(verbose){
            printf("is inf: %d\n", is_inf_gr);
            if(is_inf_gr){
//                print_cluster(cluster, cluster_size, box_size, 'd');
            }
        }
    }

    if(rng) gsl_rng_free(rng);

    return (double)res / N_iter;
}

int is_infinite_grid_v(std::vector<int>& grid, int N, std::vector<int>& is_checked, std::vector<int>& cluster, int& cluster_size)
{
    int N2 = N * N;
    int i;

    for(i = 0; i < N2; ++i) is_checked[i] = 0;

    i = 0;
    while(i < N){   // we need to check only the top row to find infinite clusters
        if(grid[i]){
            clear_cluster_v(cluster, cluster_size);
            add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, i);
            if(is_infinite_cluster_v(cluster, cluster_size, N)) {
                return 1;
            }
        } else {
            is_checked[i] = 1;
//            print_grid(is_checked, N, 'r');
        }
        ++i;
    }

    return 0;
}

int add_to_cluster_v(std::vector<int>& grid, int N, std::vector<int>& is_checked, std::vector<int>& cluster, int& cluster_size, int pos)
{
    if(!is_checked[pos]){
        is_checked[pos] = 1;
        if(grid[pos]){
            int N2 = N * N;
            cluster[cluster_size] = pos;
            ++(cluster_size);

            if(pos - N >= 0)               { add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos - N); }
            if((pos - 1) / N == (pos / N)) { add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos - 1); }
            if((pos + 1) / N == (pos / N)) { add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos + 1); }
            if(pos + N < N2)               { add_to_cluster_v(grid, N, is_checked, cluster, cluster_size, pos + N); }
        }
    }

    return 0;
}

int is_infinite_cluster_v(std::vector<int>& cluster, int& cluster_size, int N)
{
    if(cluster[0] / N > 0)
        return 0;

    int i = 0;
    int next_row_to_add = 0;
    while(i < (cluster_size)){
        if(cluster[i] / N == next_row_to_add){
            ++next_row_to_add;
            if(next_row_to_add == N)
                return 1;
        }
        ++i;
    }

    return 0;
}


int print_cluster_v(std::vector<int>& cluster, int cluster_size, int N, char prefix)
{
    int i;
    for(i = 0; i < cluster_size; ++i){
        printf("%c, %d: %d = (%d, %d)\n", prefix, i, cluster[i], cluster[i] / N, cluster[i] % N);
    }

    return 0;
}

int fill_grid_v(std::vector<int>& grid, int N, double p, gsl_rng *rng)
{
    int N2 = N * N;
    int i;

    for(i = 0; i < N2; ++i){
        grid[i] = (gsl_rng_uniform(rng) < p);
    }

    return 0;
}

int print_grid_v(std::vector<int>& grid, int N, char prefix)
{
    if(prefix > 0)
        printf("%c\n", prefix);

    int i, j;
    for(i = 0; i < N; ++i){
        for(j = 0; j < N; ++j){
            printf("%d ", grid[i * N + j]);
        }
        printf("\n");
    }

    return 0;
}

int clear_cluster_v(std::vector<int>& cluster, int& cluster_size)
{
    int i;
    for(i = 0; i < cluster_size; ++i)
        cluster[i] = -1;

    cluster_size = 0;

    return 0;
}
