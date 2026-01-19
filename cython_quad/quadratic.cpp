#include "quadratic.hpp"

int quadratic(double a, double b, double c, Roots *roots)
{
	double d = b * b - 4 * a * c;
	if(d >= 0){
        double s = sqrt(d);
        double aa = a * 2;
        roots->x1 = (s-b) / aa;
        roots->x2 = (-s-b) / aa;
        return 0;
	} else {
		printf("No real roots\n");
		return 1;
	}
}
