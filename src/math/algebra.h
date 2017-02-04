#ifndef MATH_ALGEBRA_H
#define MATH_ALGEBRA_H

//

#include "../core/param.h"

namespace math_algebra {

inline int sum_array_integer(int *array_in, int size_array_in) {

  int sum = 0;
  for (int i = 0; i < size_array_in; i++) sum += array_in[i];

  return sum;

} // sum_array_integer



inline int product_array_integer(int *array_in, int size_array_in) {

  int product = 1;
  for (int i = 0; i < size_array_in; i++) product *= array_in[i];

  return product;

} // product_array_integer



//inline int minval(int *array_in, int size_array_in) {
//
//  int minval = array_in[0];
//  for (int i = 1; i < size_array_in; i++)
//    minval = std::min(minval,array_in[i]);
//
//  return minval;
//
//} // minval
template <class T>
T minval(T *&array_in, int size_array_in) {

  T minval = array_in[0];
  for (int i = 1; i < size_array_in; i++)
    minval = std::min(minval,array_in[i]);

  return minval;

} // minval

} // math_algebra

//

#endif
