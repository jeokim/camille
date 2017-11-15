#ifndef MATH_INTERPOLATE_H
#define MATH_INTERPOLATE_H

//

#include "../core/param.h"

namespace math_interpolate {

inline double interpolate_Lagrange_1D(double *x, double *y, int num_data, double query) {

    double interpolated = 0.0;
 
    for (int i = 0; i < num_data; i++) {

      double term = y[i];
      for (int j = 0; j < num_data; j++)
        if (j != i)
          term = term * (query - x[j]) / (x[i] - x[j]);

      interpolated += term;
    } // i
 
    return interpolated;

} // interpolate_Lagrange_1D

} // math_interpolate

//

#endif
