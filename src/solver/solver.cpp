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

  } // itime_step

  return;

} // solve

} // solver
