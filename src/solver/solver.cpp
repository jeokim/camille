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
      io::write_solution_on_the_fly(myinput, mygrid, block, mystate, temporal::num_time_steps);

    } // temporal::time_step%(myinput->report_freq)
// hack for writing time-resolved pointwise data
if (temporal::time_step%(myinput->report_freq) == 0) {
std::ofstream ofs;
int l0;
double pPrime, uPrime, wplus, wminus, ws, wZ;
//
l0 = mygrid->idx1D(21-1, 0, 0);
pPrime = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]);
uPrime = mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
wplus = pPrime + uPrime;
wminus = pPrime - uPrime;
ws = mystate->sol[IVAR_S][l0]/mystate->sol_aux[IAUX_CP][l0];
wZ = mystate->sol[IVAR_Z][l0];
if (mpi::irank == 0) { 
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
ofs.open("inlet_wZ.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wZ
    << std::endl;
ofs.close();
} // mpi::irank
//
l0 = mygrid->idx1D(mygrid->ie[XI]-20, 0, 0);
pPrime = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]);
uPrime = mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
wplus = pPrime + uPrime;
wminus = pPrime - uPrime;
ws = mystate->sol[IVAR_S][l0]/mystate->sol_aux[IAUX_CP][l0];
wZ = mystate->sol[IVAR_Z][l0];
if (mpi::irank == mpi::nprocs-1) { 
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
ofs.open("outlet_wZ.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wZ
    << std::endl;
ofs.close();
} // mpi::irank
//
l0 = mygrid->idx1D(mygrid->ie[XI]-20, 4, 0);
pPrime = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]);
uPrime = mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
wplus = pPrime + uPrime;
wminus = pPrime - uPrime;
ws = mystate->sol[IVAR_S][l0]/mystate->sol_aux[IAUX_CP][l0];
wZ = mystate->sol[IVAR_Z][l0];
if (mpi::irank == mpi::nprocs-1) { 
ofs.open("outlet_r0.04_wplus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wplus
    << std::endl;
ofs.close();
ofs.open("outlet_r0.04_wminus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wminus
    << std::endl;
ofs.close();
ofs.open("outlet_r0.04_ws.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << ws
    << std::endl;
ofs.close();
ofs.open("outlet_r0.04_wZ.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wZ
    << std::endl;
ofs.close();
} // mpi::irank
//
l0 = mygrid->idx1D(mygrid->ie[XI]-20, 8, 0);
pPrime = mystate->sol[IVAR_P][l0]/(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]);
uPrime = mystate->sol[IVAR_UX][l0]/sqrt(myinput->gamma_specificheat*mystate->sol_mean[IVAR_P][l0]/mystate->sol_aux[IAUX_RHO_MEAN][l0]);
wplus = pPrime + uPrime;
wminus = pPrime - uPrime;
ws = mystate->sol[IVAR_S][l0]/mystate->sol_aux[IAUX_CP][l0];
wZ = mystate->sol[IVAR_Z][l0];
if (mpi::irank == mpi::nprocs-1) { 
ofs.open("outlet_r0.08_wplus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wplus
    << std::endl;
ofs.close();
ofs.open("outlet_r0.08_wminus.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wminus
    << std::endl;
ofs.close();
ofs.open("outlet_r0.08_ws.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << ws
    << std::endl;
ofs.close();
ofs.open("outlet_r0.08_wZ.dat", std::ofstream::app);
ofs << std::setw(16) << mygrid->cell[l0].xyz[XDIR] << " "
    << std::setw(16) << temporal::time_sol << " "
    << std::setw(16) << mygrid->cell[l0].xyz[RDIR] << " "
    << std::setw(16) << wZ
    << std::endl;
ofs.close();
} // mpi::irank
// no more hack from now on
} // itime_step%2

  } // itime_step

  return;

} // solve

} // solver
