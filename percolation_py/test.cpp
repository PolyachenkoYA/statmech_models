#include <stdio.h>
#include <stdlib.h>

#include "percol.hpp"

using namespace std;

int main(int argc, char** argv)
{
    if(argc != 8){
        printf("usage:\n%s   Ni   p1   p2   Np   BoxL   seed   id\n", argv[0]);
        return 1;
    }

    int N_iter = atoi(argv[1]);
    double p1 = atof(argv[2]);
    double p2 = atof(argv[3]);
    int Np = atoi(argv[4]);
    int box_size = atoi(argv[5]);
    int my_seed = atoi(argv[6]);
    int id = atoi(argv[7]);

    char filename[80];
    sprintf(filename, "data/Ni%d_Box%d_Seed%d_P%1.3lf_%1.3lf_%d_ID%d.dat", N_iter, box_size, my_seed, p1, p2, Np, id);

    double* p = (double*) malloc(sizeof(double) * Np);
    double* P = (double*) malloc(sizeof(double) * Np);
    for(int i = 0; i < Np; ++i) {
        p[i] = p1 + ((p2 - p1) * i) / Np;
        P[i] = prob_for_p(p[i], box_size, N_iter, my_seed);
        printf("%lf \n", (double)(i + 1) / Np);
    }

    FILE *output_file;
    output_file = fopen(filename, "w");
    for(int i = 0; i < Np; ++i) {
        fprintf(output_file, "%lf %lf\n", p[i], P[i]);
    }
    fclose(output_file);

    free(p);
    free(P);
    return 0;
}

