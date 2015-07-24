#include <math.h>
#include "cubicsolver.h"
#include "RF.h"
#define ERRORNOTCUBIC 501

//Constants for the cubic equation solver
// the function solve the equation ax^3 + bx^2 + cx + d = 0. It returns the smallest positive real root or real root
// returns 0 if the equation is trully cubic.
int cubicsolver(double a, double b, double c, double d, double roots[][2]) {

    //the equation is quadratic, but not cubic. The first coefficient must be non 0.
    if (a == 0)
    {
        SERR1("a=%e NOT OK", a);
    }


    b /= a;
    c /= a;
    d /= a;
    double disc, q, r, dum1, s, t, term1, r13;

    q = (3.0*c - (b*b))/9.0;
    r = -(27.0*d) + b*(9.0*c - 2.0*(b*b));
    r /= 54.0;
    disc = q*q*q + r*r;
    roots[0][1] = 0; //The first root is always real.
    term1 = (b/3.0);
    if (disc > 0) { // one root real, two are complex
        s = r + sqrt(disc);
        s = ((s < 0) ? -pow(-s, (1.0/3.0)) : pow(s, (1.0/3.0)));
        t = r - sqrt(disc);
        t = ((t < 0) ? -pow(-t, (1.0/3.0)) : pow(t, (1.0/3.0)));
        roots[0][0] = -term1 + s + t;
        term1 += (s + t)/2.0;
        roots[2][0] = -term1;
        roots[1][0] = -term1;
        term1 = sqrt(3.0)*(-t + s)/2;
        roots[1][1] = term1;
        roots[2][1] = -term1;
        return NOERROR;
    }
    // End if (disc > 0)
    // The remaining options are all real
    roots[1][1] = 0;
    roots[2][1] = 0;
    if (disc == 0){ // All roots real, at least two are equal.
        r13 = ((r < 0) ? -pow(-r,(1.0/3.0)) : pow(r,(1.0/3.0)));
        roots[0][0] = -term1 + 2.0*r13;
        roots[1][0] = -(r13 + term1);
        roots[2][0] = -(r13 + term1);
        return 0;
    } // End if (disc == 0)
    // Only option left is that all roots are real and unequal (to get here, q < 0)
    q = -q;
    dum1 = q*q*q;
    dum1 = acos(r/sqrt(dum1));
    r13 = 2.0*sqrt(q);
    roots[0][0]= -term1 + r13*cos(dum1/3.0);
    roots[1][0] = -term1 + r13*cos((dum1 + 2.0*M_PI)/3.0);
    roots[2][0] = -term1 + r13*cos((dum1 + 4.0*M_PI)/3.0);
    return NOERROR;
}
