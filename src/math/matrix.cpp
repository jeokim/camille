#include "matrix.h"

namespace math_matrix {

double compute_determinant(double *A, int size_of_A) {

  double determinant = 0.0;

  switch ( size_of_A ) {
  case 2:

    determinant = compute_determinant_2by2(A);

    break;

  default:

    mpi::graceful_exit("The routine computing a matrix determinant failed since computing determinant for this size of matrix is not implemented.");

    break;

  } // size_of_A

  return determinant;

} // compute_determinant



double compute_determinant_2by2(double *A) {

  return A[ij2by2[0][0]] * A[ij2by2[1][1]] - A[ij2by2[0][1]] * A[ij2by2[1][0]];

} // compute_determinant_2by2



double solve_linearSystem(double *A, double *&x, double *b, int size_of_A) {

  double determinant = 0.0;

  switch ( size_of_A ) {
  case 2:

    determinant = solve_linearSystem_2by2(A, x, b, size_of_A);

    break;

  default:

    mpi::graceful_exit("The routine solving a linear system failed since solving for this size of matrix is not implemented.");

    break;

  } // size_of_A

  return determinant;

} // solve_linearSystem



double solve_linearSystem_2by2(double *A, double *&x, double *b, int size_of_A) {

  // solve a linear system with the coefficient matrix of size 2 by 2
  // since the exact solution is known, hard-code the solution

  double determinant = compute_determinant(A, size_of_A);
  double inverse_determinant = 1.0 / determinant;

  x[0] =  A[ij2by2[1][1]] * b[0] - A[ij2by2[0][1]] * b[1];
  x[1] = -A[ij2by2[1][0]] * b[0] + A[ij2by2[0][0]] * b[1];

  for (int i = 0; i < size_of_A; i++)
    x[i] *= inverse_determinant;

  return determinant;

} // solve_linearSystem_2by2

} // math_matrix
