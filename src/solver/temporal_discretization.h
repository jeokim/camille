#ifndef SOLVER_TEMPORAL_DISCRETIZATION_H
#define SOLVER_TEMPORAL_DISCRETIZATION_H

//

#include <string>

#include "../core/param.h"
#include "../core/input.h"
#include "../math/matrix.h"
#include "../parallel/parallel.h"
#include "../core/state.h"
#include "spatial_discretization.h"
#include "../io/io.h"

namespace temporal {

class TemporalDiscretization {

  public:

    std::string name_scheme;
    int order_of_accuracy;

}; // TemporalDiscretization



class RungeKutta : public TemporalDiscretization {

  public:

    // notation following the Butcher tableau for y_{n+1} = y_n + \sum_{i = 1}^{s} b_i k_i
    //
    //   0   |
    //   c_2 | a_{21}
    //   c_3 | a_{31} a_{32}
    //   ...
    //   c_s | a_{s1} a_{s2} ... a_{s s-1}
    //   _________________________________
    //       | b_1    b_2    ... b_s
    //
    //  Thus, for a s-th order Runge--Kutta scheme, there are "s" elements for b_i; 
    //                                                        "s-1" elements for c_i; 
    //                                                        "n(n-1)/2" elements for a_ij

    int num_stages;
    double cfl_max; // maximum allowed value of CFL for a stable time integration

    double *b_i;
    double *c_i;
    double **a_ij;

    void initialize(UserInput *);
    void allocateButcherTableau(int);

    void integrate(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *, double, double, int, int, double **, double **);

}; // RungeKutta



extern std::string temporal_scheme;
extern std::string fix_dt_or_cfl;
extern double dt;
extern double *cfl;
extern double cfl_max;
//
extern int time_step_lastrun;
extern int time_step;
extern int num_time_steps;
//
extern double time_sol_lastrun;
extern double time_sol;

extern int num_vars_sol;
extern int num_samples;
extern double **rhs;
extern double **sol_interim;

extern double **yn; // solution vector at t = t_n
extern double **ynp1; // solution vector at t = t_(n+1)

extern RungeKutta rk;



void initialize(UserInput *, State *);
void finalize(void);

void precompute_something(UserInput *, Geometry::StructuredGrid *, State *);
void time_advance(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *); // advance a solution at a time leven n to n+1

void iblank_correction(int, int, double **, Geometry::StructuredGrid *, double []);

double compute_CFL(UserInput *, Geometry::StructuredGrid *, State *);

} // temporal

//

#endif
