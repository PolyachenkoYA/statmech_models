#include <stdio.h>
#include <math.h>

typedef struct {
    double x1, x2;
} Roots;

int quadratic(double a, double b, double c, Roots *roots);
