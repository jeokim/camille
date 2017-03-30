#include <iostream>
#include <iomanip>
#include <fstream>

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

    } // temporal::time_step%(myinput->report_freq)
    io::write_solution_on_the_fly(myinput, mygrid, block, mystate, temporal::num_time_steps);

// hack for writing time-resolved pointwise data
if (temporal::time_step%(myinput->report_freq) == 0) {
std::ofstream ofs;
int i, l0;
double pPrime, uPrime, wplus, wminus, ws, wZ;
//
i = NONE;
for (int i_query = mygrid->is[XI]; i_query <= mygrid->ie[XI]; i_query++) {
  l0 = mygrid->idx1D(i_query, mygrid->is[ETA], mygrid->is[ZETA]);
  if (mygrid->cell[l0].xyz[XDIR] == 0.0) // at the inlet
    i = i_query;
} // i_query
if (i != NONE) {
wplus = 0.0; wminus = 0.0; ws = 0.0; wZ = 0.0;
for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
l0 = mygrid->idx1D(i, j, 0);
pPrime = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]);
uPrime = mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
wplus += pPrime + uPrime;
wminus += pPrime - uPrime;
ws += mystate->sol[IVAR_S][l0]/(mystate->sol_aux[IAUX_CP])[l0];
if (myinput->model_pde == "LEE_MIXFRAC_CONSTGAMMA")
wZ += mystate->sol[IVAR_Z][l0];
} // j
wplus /= mygrid->num_cells_dir[ETA]; wminus /= mygrid->num_cells_dir[ETA]; ws /= mygrid->num_cells_dir[ETA]; wZ /= mygrid->num_cells_dir[ETA];
if (temporal::time_sol >= 4.0) { // after transients
ofs.open("inlet_wplus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wplus
    << std::endl;
ofs.close();
ofs.open("inlet_wminus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wminus
    << std::endl;
ofs.close();
ofs.open("inlet_ws.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << ws
    << std::endl;
ofs.close();
if (myinput->model_pde == "LEE_MIXFRAC_CONSTGAMMA") {
ofs.open("inlet_wZ.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wZ
    << std::endl;
ofs.close();
} // myinput->model_pde
} // temporal::time_sol
} // i
//
i = NONE;
for (int i_query = mygrid->is[XI]; i_query <= mygrid->ie[XI]; i_query++) {
  l0 = mygrid->idx1D(i_query, mygrid->is[ETA], mygrid->is[ZETA]);
  if (mygrid->cell[l0].xyz[XDIR] == 1.0) // at the outlet
    i = i_query;
} // i_query
if (i != NONE) {
wplus = 0.0; wminus = 0.0; ws = 0.0; wZ = 0.0;
for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
l0 = mygrid->idx1D(i, j, 0);
pPrime = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]);
uPrime = mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
wplus += pPrime + uPrime;
wminus += pPrime - uPrime;
ws += mystate->sol[IVAR_S][l0]/(mystate->sol_aux[IAUX_CP])[l0];
if (myinput->model_pde == "LEE_MIXFRAC_CONSTGAMMA")
wZ += mystate->sol[IVAR_Z][l0];
} // j
wplus /= mygrid->num_cells_dir[ETA]; wminus /= mygrid->num_cells_dir[ETA]; ws /= mygrid->num_cells_dir[ETA]; wZ /= mygrid->num_cells_dir[ETA];
if (temporal::time_sol >= 4.0) { // after transients
ofs.open("outlet_wplus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wplus
    << std::endl;
ofs.close();
ofs.open("outlet_wminus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wminus
    << std::endl;
ofs.close();
ofs.open("outlet_ws.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << ws
    << std::endl;
ofs.close();
if (myinput->model_pde == "LEE_MIXFRAC_CONSTGAMMA") {
ofs.open("outlet_wZ.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wZ
    << std::endl;
ofs.close();
} // myinput->model_pde
} // temporal::time_sol
} // i
// no more hack from now on
} // temporal::time_step%(myinput->report_freq)

  } // itime_step

  return;

} // solve

} // solver
