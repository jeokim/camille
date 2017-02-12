#ifndef CORE_MACROSINLINES_H
#define CORE_MACROSINLINES_H

//

#include <string.h>
#include <string>
#include <algorithm>
#include <functional> 
#include <cctype>
#include <locale>
#include "param.h"
#include <assert.h>

// integer array
#define ALLOCATE1D_INT_1ARG(PTR1D,N1) {PTR1D = new int[N1]; for (int _i = 0; _i < N1; _i++) PTR1D[_i] = 0;}
#define ALLOCATE1D_INT_2ARG(PTR1D,N1,N2) {PTR1D = new int[N1 * N2]; for (int _i = 0; _i < N1 * N2; _i++) PTR1D[_i] = 0;}
#define ALLOCATE2D_INT(PTR2D,N1,N2) {PTR2D = new int*[N1]; for (int _i = 0; _i < N1; _i++) {PTR2D[_i] = new int[N2]; for (int _j = 0; _j < N2; _j++) PTR2D[_i][_j] = 0;}}
//
// double array
#define ALLOCATE1D_DOUBLE_1ARG(PTR1D,N1) {PTR1D = new double[N1]; for (int _i = 0; _i < N1; _i++) PTR1D[_i] = 0.0;}
#define ALLOCATE1D_DOUBLE_2ARG(PTR1D,N1,N2) {PTR1D = new double[N1 * N2]; for (int _i = 0; _i < N1 * N2; _i++) PTR1D[_i] = 0.0;}
#define ALLOCATE2D_DOUBLE(PTR2D,N1,N2) {PTR2D = new double*[N1]; for (int _i = 0; _i < N1; _i++) {PTR2D[_i] = new double[N2]; for (int _j = 0; _j < N2; _j++) PTR2D[_i][_j] = 0.0;}}
//
// deallocate arrays
#define DEALLOCATE_1DPTR(PTR1D) {delete[] PTR1D;}
#define DEALLOCATE_2DPTR(PTR2D,N1) {for (int _i = 0; _i < N1; _i++) delete[] PTR2D[_i]; DEALLOCATE_1DPTR(PTR2D);}
#define DEALLOCATE_3DPTR(PTR3D,N1,N2) {for (int _i = 0; _i < N1; _i++) for (int _j = 0; _j < N2; _j++) delete[] PTR3D[_i][_j]; DEALLOCATE_2DPTR(PTR3D,N1);}
#define DEALLOCATE_4DPTR(PTR4D,N1,N2,N3) {for (int _i = 0; _i < N1; _i++) for (int _j = 0; _j < N2; _j++) for (int _k = 0; _k < N3; _k++) delete[] PTR4D[_i][_j][_k]; DEALLOCATE_3DPTR(PTR4D,N1,N2);}
#define DEALLOCATE_5DPTR(PTR5D,N1,N2,N3,N4) {for (int _i = 0; _i < N1; _i++) for (int _j = 0; _j < N2; _j++) for (int _k = 0; _k < N3; _k++) for (int _l = 0; _l < N4; _l++) delete[] PTR5D[_i][_j][_k][_l]; DEALLOCATE_4DPTR(PTR5D,N1,N2,N3);}

#define MESSAGE_STDOUT(M) {if (mpi::irank == 0) std::cout << ">>> " << M << std::endl;}



namespace todos {

extern int num_todos;

void add(std::string);

} // todos



inline const char *cstr_to_constchar(std::string cstr) {

  const char *constchar = cstr.c_str();

  return constchar;

} // c str_to_constchar



inline int maxloc(int *val, int valsize) {

  int imax = 0;
  int max = val[imax];

  for (int i = 1; i < valsize; i++) {

    if (val[i] > max) {

      imax = i;
      max = val[i];

    }
  }

  return imax;

} // maxloc



inline void str_tolower(std::string &str_in_n_out) {

  for (int i = 0; i < str_in_n_out.size(); i++)
    str_in_n_out[i] = std::tolower(str_in_n_out[i]);

  return;

} // str_tolower



inline void str_toupper(std::string &str_in_n_out) {

  for (int i = 0; i < str_in_n_out.size(); i++)
    str_in_n_out[i] = std::toupper(str_in_n_out[i]);

  return;

} // str_toupper



// following trimming functions are written by Evan Teran (http://stackoverflow.com/)
// trim from start (in place)
static inline void LTRIM(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void RTRIM(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void TRIM(std::string &s) {
    LTRIM(s);
    RTRIM(s);
}

// trim from start (copying)
static inline std::string LTRIMMED(std::string s) {
    LTRIM(s);
    return s;
}

// trim from end (copying)
static inline std::string RTRIMMED(std::string s) {
    RTRIM(s);
    return s;
}

// trim from both ends (copying)
static inline std::string TRIMMED(std::string s) {
    TRIM(s);
    return s;
}
// Thanks, Evan Teran! (http://stackoverflow.com/)

//

#endif
