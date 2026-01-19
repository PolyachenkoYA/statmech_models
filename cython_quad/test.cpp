#include "quadratic.hpp"

int main()
{
    Roots r;
    double a,b,c;

    printf("enter a,b,c:\n");
    scanf("%lf%lf%lf", &a, &b, &c);

    quadratic(a, b, c, &r);

    printf("x1,2: %lf, %lf\n", r.x1, r.x2);

    return 0;
}
