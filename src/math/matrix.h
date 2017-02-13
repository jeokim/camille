#ifndef MATH_MATRIX_H
#define MATH_MATRIX_H

//

#include <cmath>

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"

namespace math_matrix {

// determinant
double compute_determinant(double *, int);
double compute_determinant_2by2(double *);

// solver
double solve_linearSystem(double *, double *&, double *, int);
double solve_linearSystem_2by2(double *, double *&, double *, int);



inline double inner_product(double *vector0, double *vector1, int size_vector) {

  double product = 0.0;

  for (int i = 0; i < size_vector; i++)
    product += vector0[i] * vector1[i];

  return product;

} // inner_product



inline double cross_product(double *vector0, double *vector1, int size_vector) {

  if (size_vector != 2)
    mpi::graceful_exit("A cross product is attempted for a vector of unknown size.");

  return vector0[0] * vector1[1] - vector0[1] * vector1[0];

} // cross_product



inline double cross_product_normalized(double *vector0, double *vector1, int size_vector) {

  if (size_vector != 2)
    mpi::graceful_exit("A cross product is attempted for a vector of unknown size.");

  double product0 = 0.0, product1 = 0.0;

  for (int i = 0; i < size_vector; i++) {
    product0 += vector0[i] * vector0[i];
    product1 += vector1[i] * vector1[i];
  } // i

  return (vector0[0] * vector1[1] - vector0[1] * vector1[0]) / sqrt(product0 * product1);

} // cross_product

} // math_matrix

//

#endif
