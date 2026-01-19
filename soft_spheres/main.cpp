#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "soft_spheres.h"

int main(int argc, char **argv) {
    if(argc != 7){
        printf("usage:\n%s   L/a   n   Nt   d   Temp   seed\n", argv[0]);
        return  1;
    }

//    int s = 2; // L = s * a
//    double n = 0.5;
//    int Nt = 100;
//    int my_seed = 0;
    int s = atoi(argv[1]);
    double n = atof(argv[2]);
    int Nt = atoi(argv[3]);
    double d = atof(argv[4]);
    double Temp = atof(argv[5]);
    int my_seed = atoi(argv[6]);

    double* U = (double *) malloc(sizeof (double) * Nt);
    double* P = (double *) malloc(sizeof (double) * Nt);

    compute_E_evol_C(s, n, Nt, d, Temp, my_seed, U, P, 0);

    char filename[180];
    sprintf(filename, "s%d_n%lf_Nt%d_d%lf_Temp%lf_seed%d.dat", s, n, Nt, d, Temp, my_seed);

    FILE *output_file;
    output_file = fopen(filename, "w");
    for(int it = 0; it < Nt; ++it) {
//        fprintf(output_file, "%lf %lf\n", U[it], U2[it]);
        fprintf(output_file, "%lf %lf\n", U[it], P[it]);
    }
    fclose(output_file);
    printf("saved to '%s'\n", filename);

    free(U);
    free(P);
    printf("freed\n");

    return 0;
}
