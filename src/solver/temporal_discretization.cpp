#include <cmath>
#include <iostream>

#include "temporal_discretization.h"

namespace temporal {

std::string temporal_scheme;
std::string fix_dt_or_cfl;
double dt;
double *cfl;
double cfl_max;
//
int time_step_lastrun;
int time_step;
int num_time_steps;
//
double time_sol_lastrun;
double time_sol;

int num_vars_sol;
int num_samples;
double **rhs;
double **sol_interim;

double **yn; // solution vector at t = t_n
double **ynp1; // solution vector at t = t_(n+1)

RungeKutta rk;



void RungeKutta::initialize(UserInput *myinput) {

  if (temporal_scheme == "RUNGE_KUTTA2") {
    name_scheme = "2nd-order standard Runge--Kutta";
    order_of_accuracy = 2;

    num_stages = 2;
    //cfl_max = ... is what?

    // build the Butcher tableau
    allocateButcherTableau(order_of_accuracy);
    //
    b_i[FIRST] = 1.0 / 2.0;
    b_i[SECOND] = 1.0 / 2.0;
    //
    c_i[FIRST] = 0.0;
    c_i[SECOND] = 1.0;
    //
    a_ij[SECOND][FIRST] = 1.0;

  } // temporal_scheme
  else if (temporal_scheme == "RUNGE_KUTTA4") {

    name_scheme = "4th-order standard Runge--Kutta";
    order_of_accuracy = 4;

    num_stages = 4;
    //cfl_max = ... is what?

    // build the Butcher tableau
    allocateButcherTableau(order_of_accuracy);
    //
    b_i[FIRST] = 1.0 / 6.0;
    b_i[SECOND] = 1.0 / 3.0;
    b_i[THIRD] = 1.0 / 3.0;
    b_i[FOURTH] = 1.0 / 6.0;
    //
    c_i[FIRST] = 0.0;
    c_i[SECOND] = 1.0 / 2.0;
    c_i[THIRD] = 1.0 / 2.0;
    c_i[FOURTH] = 1.0;
    //
    a_ij[SECOND][FIRST] = 1.0 / 2.0;
    a_ij[THIRD][SECOND] = 1.0 / 2.0;
    a_ij[FOURTH][THIRD] = 1.0;

  } // temporal_scheme
  else
    mpi::graceful_exit("Unknown temporal scheme.");

  return;

} // RungeKutta::initialize



void RungeKutta::allocateButcherTableau(int order_of_accuracy_in) {

  ALLOCATE1D_DOUBLE_1ARG(b_i, order_of_accuracy_in);
  ALLOCATE1D_DOUBLE_1ARG(c_i, order_of_accuracy_in);
  ALLOCATE2D_DOUBLE(a_ij, order_of_accuracy_in, order_of_accuracy_in);

  for (int i = FIRST; i < order_of_accuracy_in; i++)
    b_i[i] = 0.0;

  for (int i = FIRST; i < order_of_accuracy_in; i++)
    c_i[i] = 0.0;

  for (int i = FIRST; i < order_of_accuracy_in; i++)
    for (int j = FIRST; j < order_of_accuracy_in; j++)
      a_ij[i][j] = 0.0;

  return;

} // RungeKutta::initialize



//void RungeKutta::integrate(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate, 
//                           double tn, double h, int num_vars, int num_samples, double *yn[], double *ynp1[]) {
//
//  // this solves an initial value problem, y'(t) = f(t, y(t)) with y(0) = y_0, discretized as 
//  //   y_(n+1) = y_n + h (k_1/2 + k_2/2)
//  // where
//  //   k_1 = f(t_n,     y_n)
//  //   k_2 = f(t_n + h, y_n + h k_1)
//
//  // y_(n+1) = y_n
//  for (int ivar = 0; ivar < num_vars; ivar++)
//    for (int l0 = 0; l0 < num_samples; l0++)
//      (ynp1[ivar])[l0] = (yn[ivar])[l0];
//
//  // 1st predicting stage
//  spatial::compute_RHS(myinput, mygrid, mystate, yn, rhs, tn + c_i[FIRST] * h); // k_1
//  //
//  for (int ivar = 0; ivar < num_vars; ivar++)
//    for (int l0 = 0; l0 < num_samples; l0++) {
//
//      if (mygrid->cell[l0].iblank != BLANKED) { // only non-hole cells are time advanced
//
//        (rhs[ivar])[l0] *= h; // h k_1
//        (ynp1[ivar])[l0] += b_i[FIRST] * (rhs[ivar])[l0]; // h (k_1 / 2)
//
//        (sol_interim[ivar])[l0] = (yn[ivar])[l0] + a_ij[SECOND][FIRST] * (rhs[ivar])[l0]; // y_n + h k_1
//
//      } // mygrid->cell[l0].iblank
//      else { // hole points bypass
//
//        (rhs[ivar])[l0] = 0.0;
//        (ynp1[ivar])[l0] = mystate->sol_ref[ivar];
//
//        (sol_interim[ivar])[l0] = mystate->sol_ref[ivar];
//
//      } // mygrid->cell[l0].iblank
//    } // l0
//  spatial::update_boundary(myinput, mygrid, mystate, num_vars, sol_interim, tn + c_i[SECOND] * h); // apply boundary conditions for the next stage
//
//  // 2nd Runge--Kutta correcting stage
//  spatial::compute_RHS(myinput, mygrid, mystate, sol_interim, rhs, tn + c_i[SECOND] * h); // k_2
//  //
//  for (int ivar = 0; ivar < num_vars; ivar++)
//    for (int l0 = 0; l0 < num_samples; l0++) {
//
//      if (mygrid->cell[l0].iblank != BLANKED) { // only non-hole cells are time advanced
//
//        (rhs[ivar])[l0] *= h; // h k_2
//        (ynp1[ivar])[l0] += b_i[SECOND] * (rhs[ivar])[l0]; // h (k_2 / 2)
//
//      } // mygrid->cell[l0].iblank
//      else { // hole points bypass
//
//        (rhs[ivar])[l0] = 0.0;
//        (ynp1[ivar])[l0] = mystate->sol_ref[ivar];
//
//      } // mygrid->cell[l0].iblank
//    } // l0
//  spatial::update_boundary(myinput, mygrid, mystate, num_vars, ynp1, tn + c_i[SECOND] * h); // apply boundary conditions for the next stage
//
//  return;
//
//} // RungeKutta::integrate



void RungeKutta::integrate(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate, 
                           double tn, double h, int num_vars, int num_samples_in, double **y_cur, double **y_next) {

  // this function is tested for cases solving an initial value problem, y'(t) = f(t, y(t)) with y(0) = y_0, discretized by 
  // 2nd order standard explict Runge--Kutta 
  //   y_(n+1) = y_n + h (k_1/2 + k_2/2)
  // where
  //   k_1 = f(t_n,     y_n)
  //   k_2 = f(t_n + h, y_n + h k_1)
  // and
  // 4th order standard explict Runge--Kutta 
  //   y_(n+1) = y_n + h (k_1/6 + k_2/3 + k_3/3 + k_4/6)
  // where
  //   k_1 = f(t_n,       y_n)
  //   k_2 = f(t_n + h/2, y_n + h/2 k_1)
  //   k_3 = f(t_n + h/2, y_n + h/2 k_2)
  //   k_4 = f(t_n + h,   y_n + h   k_3)

  for (int ivar = 0; ivar < num_vars; ivar++)
    for (int l0 = 0; l0 < num_samples_in; l0++) {

      (y_next[ivar])[l0] = (y_cur[ivar])[l0]; // y_(n+1) = y_n
      (sol_interim[ivar])[l0] = (y_cur[ivar])[l0]; // initial predictor is equal to the last time step's converged solution

    } // l0

  int last_stage = num_stages - 1;
  for (int istage = FIRST; istage <= last_stage; istage++) {

    spatial::compute_RHS(myinput, mygrid, mystate, sol_interim, rhs, tn + c_i[istage] * h); // k_i
    //
    for (int l0 = 0; l0 < num_samples_in; l0++) {
      if (mygrid->cell[l0].iblank != BLANKED) { // only non-hole cells are time advanced
        for (int ivar = 0; ivar < num_vars; ivar++) {

          (rhs[ivar])[l0] *= h; // h k_i
          (y_next[ivar])[l0] += b_i[istage] * (rhs[ivar])[l0]; // h b_i k_i; corrector

          if (istage < last_stage) // if not a last stage, compute a predictor
            (sol_interim[ivar])[l0] = (y_cur[ivar])[l0] + a_ij[istage + 1][istage] * (rhs[ivar])[l0]; // y_n + a_ij h k_i

        } // ivar
      } // mygrid->cell[l0].iblank
      else { // hole points bypass
        for (int ivar = 0; ivar < num_vars; ivar++) {

          (rhs[ivar])[l0] = 0.0;
          (y_next[ivar])[l0] = mystate->sol_ref[ivar];

          (sol_interim[ivar])[l0] = mystate->sol_ref[ivar];

        } // ivar
      } // mygrid->cell[l0].iblank
    } // l0

    // apply boundary conditions for the next stage
    if (istage < last_stage) {

      spatial::update_boundary(myinput, mygrid, mystate, num_vars, sol_interim, tn + c_i[std::min(istage + 1, last_stage)] * h);
      mystate->compute_dependent_variables(sol_interim);

    } // istage
    else {

      spatial::update_boundary(myinput, mygrid, mystate, num_vars, y_next, tn + c_i[std::min(istage + 1, last_stage)] * h);
      mystate->compute_dependent_variables(y_next);

    } // istage

    mpi::wait_allothers();

  } // istage

  return;

} // RungeKutta::integrate



void initialize(UserInput *myinput, State *mystate) {

  temporal_scheme = myinput->temporal_scheme;
  fix_dt_or_cfl = myinput->fix_dt_or_cfl;
  dt = myinput->dt;
  ALLOCATE1D_DOUBLE_1ARG(cfl, mystate->num_samples);
  for (int l0 = 0; l0 < mystate->num_samples; l0++)
    cfl[l0] = myinput->cfl;
  cfl_max = myinput->cfl;
  //
  time_step_lastrun = mystate->time_step;
  time_step = time_step_lastrun;
  num_time_steps = myinput->num_time_steps;
  //
  time_sol_lastrun = mystate->time_sol;
  time_sol = time_sol_lastrun;

  num_vars_sol = mystate->num_vars_sol;
  num_samples = mystate->num_samples;
  rhs = new double *[num_vars_sol];
  sol_interim = new double *[num_vars_sol];
  for (int ivar = 0; ivar < num_vars_sol; ivar++) {

    rhs[ivar] = new double[mystate->num_samples];
    sol_interim[ivar] = new double[mystate->num_samples];

  } // ivar

  std::string name_scheme;

  if (temporal_scheme == "RUNGE_KUTTA2" || temporal_scheme == "RUNGE_KUTTA4") {

    rk.initialize(myinput);

    name_scheme = rk.name_scheme;

    yn = new double *[num_vars_sol];
    ynp1 = new double *[num_vars_sol];
    for (int ivar = 0; ivar < num_vars_sol; ivar++) {

      yn[ivar] = new double[num_samples];
      ynp1[ivar] = new double[num_samples];

    } // ivar

  } // temporal_scheme
  else
    mpi::graceful_exit("Unknown temporal scheme.");
  //
  mpi::wait_allothers("Temporal discretization is " + name_scheme + " and is initialized.");

  return;

}



void finalize() {

  DEALLOCATE_2DPTR(rhs, num_vars_sol);
  DEALLOCATE_2DPTR(sol_interim, num_vars_sol);

  return;

} // finalize



void precompute_something(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate) {



  return;

} // precompute_something



void time_advance(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate) {

//  if (myinput->do_filter == TRUE) {
//    spatial::apply_filter(num_vars_sol, num_samples, mystate->sol, mygrid);
//    return;
//  } // myinput->do_filter

  for (int ivar = 0; ivar < num_vars_sol; ivar++) {
    for (int l0 = 0; l0 < num_samples; l0++) {

      (yn[ivar])[l0] = (mystate->sol_old[ivar])[l0];
      (ynp1[ivar])[l0] = DUMMY_DOUBLE;

    } // l0
  } // ivar

  if (temporal_scheme == "RUNGE_KUTTA2" || temporal_scheme == "RUNGE_KUTTA4")
    rk.integrate(myinput, mygrid, block, mystate, time_sol, dt, num_vars_sol, num_samples, yn, ynp1);

  else
    mpi::graceful_exit("Unknown temporal scheme.");

  // apply a selective filter
  if (myinput->do_filter == TRUE)
    spatial::apply_filter(num_vars_sol, num_samples, ynp1, mygrid);

  // correct based upon grid-blank (mask) information
  iblank_correction(num_vars_sol, num_samples, ynp1, mygrid, mystate->sol_ref);

  // update the current solution
  for (int ivar = 0; ivar < num_vars_sol; ivar++) {
    for (int l0 = 0; l0 < num_samples; l0++) {

      (mystate->sol[ivar])[l0] = (ynp1[ivar])[l0];

    } // l0
  } // ivar

  mpi::wait_allothers();

  return;

} // time_advance



void iblank_correction(int num_vars_in, int num_samples_in, double **y, Geometry::StructuredGrid *mygrid, double sol_ref[]) {

  for (int ivar = 0; ivar < num_vars_in; ivar++) {
    for (int l0 = 0; l0 < num_samples_in; l0++) {

      if (mygrid->cell[l0].iblank == BLANKED)
        (y[ivar])[l0] = sol_ref[ivar]; // if a cell is blanked out (masked), its solutions are fixed as those of a prescribed ambient state

    } // l0
  } // ivar

  return;

} // iblank_correction



double compute_CFL(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate) {

  double cfl_max_local = -DUMMY_LARGE;
  double cfl_max_global = -DUMMY_LARGE;

  double velocity_Cartesian[DIM_MAX];
  double velocity_contravariant[DIM_MAX];
  double gammaMinus1 = mystate->gamma_specificheat - 1.0;

  if (myinput->model_pde == "LINEAR_EULER") {

    for (int l0 = 0; l0 < num_samples; l0++) {

      cfl[l0] = 0.0;

      if (mygrid->cell[l0].iblank == NOTBLANKED) {

        // get a contravariant velocity
        for (int idir = XDIR; idir < DIM_MAX; idir++)
          velocity_Cartesian[idir] = (mystate->sol_mean[IVAR_UX + idir])[l0]; //+ (mystate->sol[IVAR_UX + idir])[l0];
        for (int idir = myinput->num_dim; idir < DIM_MAX; idir++)
          velocity_Cartesian[idir] = 0.0; // just in case
        mystate->to_contravariant_velocity(mygrid->cell[l0].metrics, velocity_Cartesian, velocity_contravariant);

        // contribution from flow convection
        for (int icontra = XI; icontra < myinput->num_dim; icontra++)
          for (int iCart = XDIR; iCart < myinput->num_dim; iCart++)
            cfl[l0] += abs(velocity_contravariant[icontra] * mygrid->cell[l0].metrics[icontra][iCart]);

        // contribution from acoustic propagation
        double fac = 0.0;
        for (int icontra = XI; icontra < myinput->num_dim; icontra++)
          for (int iCart = XDIR; iCart < myinput->num_dim; iCart++)
            fac +=  mygrid->cell[l0].metrics[icontra][iCart] * mygrid->cell[l0].metrics[icontra][iCart];
        double speed_of_sound = sqrt(gammaMinus1 * (mystate->sol_aux[IAUX_T_MEAN])[l0]);
        cfl[l0] += speed_of_sound * sqrt(fac);

        // cfl
        cfl[l0] *= dt * mygrid->cell[l0].Jac;
        cfl_max_local = std::max(cfl_max_local, cfl[l0]);

      } // mygrid->cell[l0].iblank
    } // l0
  } // myinput->model_pde
  else
    mpi::graceful_exit("CFL computation for other than the linearized Euler model is not implemented.");

  // get a global cfl_max
  MPI_Allreduce(&cfl_max_local, &cfl_max_global, 1, MPI_DOUBLE, MPI_MAX, mpi::comm_region);

  return cfl_max_global;

} // compute_CFL

} // temporal
