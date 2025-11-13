#ifndef POTENTIAL_H
#define POTENTIAL_H

/* 
 * Potencial de un pozo finito en 1D:
 *  V(x) = -V_0  si |x| < a
 *  V(x) =  0   si |x| >= a
 */

double V(double x, double a, double V0);

#endif

