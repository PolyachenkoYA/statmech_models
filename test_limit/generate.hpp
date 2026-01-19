#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>


typedef struct {
    double x1, x2;
} Roots;

int quadratic(double a, double b, double c, Roots *roots);
