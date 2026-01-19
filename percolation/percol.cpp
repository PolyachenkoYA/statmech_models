//
// Created by ypolyach on 10/15/21.
//

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include "percol.hpp"

double prob_for_p(double p, int box_size, int N_iter, int my_seed, bool verbose)
{
    int box_size2 = box_size * box_size;
    int* grid  = (int*)malloc(sizeof(int) * box_size2);
    int* cluster = (int*)malloc(sizeof(int) * box_size2);
    int* is_checked = (int*)malloc(sizeof(int) * box_size2);

    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, my_seed);

    int i;
    int cluster_size = box_size * box_size;   // initial N2 so it's all initialized to -1
    int res = 0;
    int is_inf_gr;

    for(i = 0; i < N_iter; ++i){
        fill_grid(grid, box_size, p, rng);

        if(verbose) {
            print_grid(grid, box_size);
        }

        is_inf_gr = is_infinite_grid(grid, box_size, is_checked, cluster, &cluster_size);
        res += is_inf_gr;

        if(verbose){
            printf("is inf: %d\n", is_inf_gr);
            if(is_inf_gr){
//                print_cluster(cluster, cluster_size, box_size, 'd');
            }
        }
    }

//    if(grid) free(grid);
//    if(is_checked) free(is_checked);
//    if(cluster) free(cluster);
//    if(rng) gsl_rng_free(rng);

    return (double)res / N_iter;
}

int is_infinite_grid(int* grid, int N, int* is_checked, int* cluster, int* cluster_size)
{
    int N2 = N * N;
    int i;

    for(i = 0; i < N2; ++i) is_checked[i] = 0;

    i = 0;
    while(i < N){   // we need to check only the top row to find infinite clusters
        if(grid[i]){
            clear_cluster(cluster, cluster_size);
            add_to_cluster(grid, N, is_checked, cluster, cluster_size, i);
            if(is_infinite_cluster(cluster, cluster_size, N)) {
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

int is_infinite_cluster(int* cluster, int* cluster_size, int N)
{
//    if(cluster[0] >= N) // we start only with 1st row so clusters will start with cluster[0] < N
//        return 0;

    int i = 0;
    int next_row_to_add = 0;
    while(i < (*cluster_size)){
        if(cluster[i] / N == next_row_to_add){
            ++next_row_to_add;
            if(next_row_to_add == N)
                return 1;
        }
        ++i;
    }

    return 0;
}

int add_to_cluster(int* grid, int N, int* is_checked, int* cluster, int* cluster_size, int pos)
{
//    return 1;
    if(!is_checked[pos]){
        is_checked[pos] = 1;
        if(grid[pos]){
            int N2 = N * N;
            cluster[*cluster_size] = pos;
            ++(*cluster_size);

            if(pos - N >= 0) add_to_cluster(grid, N, is_checked, cluster, cluster_size, pos - N);
            if((pos - 1) / N == (pos / N)) add_to_cluster(grid, N, is_checked, cluster, cluster_size, pos - 1);
            if((pos + 1) / N == (pos / N)) add_to_cluster(grid, N, is_checked, cluster, cluster_size, pos + 1);
            if(pos + N < N2) add_to_cluster(grid, N, is_checked, cluster, cluster_size, pos + N);
        }
    }

    return 0;
}

int print_cluster(int* cluster, int cluster_size, int N, char prefix)
{
    int i;
    for(i = 0; i < cluster_size; ++i){
        printf("%c, %d: %d = (%d, %d)\n", prefix, i, cluster[i], cluster[i] / N, cluster[i] % N);
    }

    return 0;
}

int fill_grid(int* grid, int N, double p, gsl_rng *rng)
{
    int N2 = N * N;
    int i;

    for(i = 0; i < N2; ++i){
        grid[i] = (gsl_rng_uniform(rng) < p);
    }

    return 0;
}

int print_grid(int* grid, int N, char prefix)
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

int clear_cluster(int* cluster, int *cluster_size)
{
    int i;
    for(i = 0; i < *cluster_size; ++i)
        cluster[i] = -1;

    *cluster_size = 0;

    return 0;
}

int copy_arr(int* from, int* to, int N)
{
    for(int i = 0; i < N; ++i) to[i] = from[i];
    return 0;
}
