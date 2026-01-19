#include "test.cpp"

int main(int argc, char** argv) {
    if(argc != 5){
        printf("usage:\n%s   N   T   Nt   seed\n");
        return 1;
    }

    int N = atoi(argv[1]);
    double Temp = atof(argv[2]);
    int Nt = atoi(argv[3]);
    int my_seed = atoi(argv[4]);
    int i;

    double *E = (double*) malloc(sizeof(double) * Nt);
    compute_evolution_C(N, Temp, Nt, my_seed, E);

    char filename[80];
    sprintf(filename, "N%d_T%lf_Nt%d_seed%d.dat", N, Temp, Nt, my_seed);

    FILE *output_file;
    output_file = fopen(filename, "w");
    for(i = 0; i < Nt; ++i) {
        fprintf(output_file, "%lf ", E[i]);
    }
    fclose(output_file);

    free(E);

    return 0;
}


