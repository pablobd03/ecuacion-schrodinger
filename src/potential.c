#include "potential.h"

/*
 * Pozo finito:
 * V(x) = -V0 dentro del pozo |x| < a
 * V(x) = 0 afuera
 */

double V(double x, double a, double V0) {
    if (x > -a && x < a)
        return -V0;
    return 0.0;
}
