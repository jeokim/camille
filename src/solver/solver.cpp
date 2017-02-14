#include <iostream>
#include <iomanip>

#include "solver.h"

namespace solver {

void initialize(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate) {

  temporal::initialize(myinput, mystate);
  spatial::initialize(myinput, mygrid, mystate);

  return;

} // initialize



void finalize() {

  temporal::finalize();
  spatial::finalize();

  return;

} // finalize



void precompute_something(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate) {

  temporal::precompute_something(myinput, mygrid, mystate);
  spatial::precompute_something(myinput, mygrid, mystate);

  return;

} // precompute_something



void solve(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate) {

  for (int itime_step = temporal::time_step_lastrun + 1; itime_step < temporal::num_time_steps + 1; itime_step++) {

    // check user input here: user inputs include which kinds of actions, which time step are they applied, and so on
    // dump the last time-step's solutions to "old" arrays
    mystate->backup_the_current_solution();

    // time-step advancement
    temporal::time_step = itime_step;
    mystate->time_step = temporal::time_step;

    // compute CFL
    double cfl_max = temporal::compute_CFL(myinput, mygrid, mystate);

    // adjust dt based upon the maximum CFL, if necessary
    double dt = temporal::dt;
    if (temporal::fix_dt_or_cfl == "FIX_CFL") {

      double ratio = temporal::cfl_max / cfl_max;
      dt *= ratio;
      temporal::dt = dt;

    } // temporal::fix_dt_or_cfl
    else if (temporal::fix_dt_or_cfl == "FIX_DT") {

      temporal::cfl_max = cfl_max;

    } // temporal::fix_dt_or_cfl

    // time advancement
    temporal::time_sol += dt;
    mystate->time_sol = temporal::time_sol;
    //
    temporal::time_advance(myinput, mygrid, block, mystate);

    // see what we have done
    if (temporal::time_step%(myinput->report_freq) == 0) {

      io::report_timeadvancing(temporal::time_step, temporal::time_sol, dt, cfl_max, myinput, mygrid, mystate);
      io::write_solution_on_the_fly(myinput, mygrid, block, mystate, temporal::num_time_steps);

    } // temporal::time_step%(myinput->report_freq)

if (itime_step%2 == 0) {
std::ofstream ofs;
//
int l0 = mygrid->idx1D(8-1, 0, 0);
double wplus = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]) + 
               mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
double ws = mystate->sol[IVAR_S][l0];
ofs.open("inlet.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR]
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR]
    << std::setw(16) << wplus
    << std::setw(16) << ws
    << std::endl;
ofs.close();
//
int l0 = mygrid->idx1D(180-1, 0, 1);
double wplus = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]) + 
               mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
double ws = mystate->sol[IVAR_S][l0];
ofs.open("outlet.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR]
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR]
    << std::setw(16) << wplus
    << std::setw(16) << ws
    << std::endl;
ofs.close();
} // itime_step%2

  } // itime_step

  return;

} // solve

} // solver
