#include "boundary_condition.h"
#include "source.h"

#include <cmath>
#include <iostream>
#include <iomanip>

namespace source {

void add_RHSSourceTerms(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs, double time) {

  apply_physicalSource(myinput, mygrid, mystate, y, rhs, time); // add physical source terms to the governing equations

  apply_bufferZone(myinput, mygrid, mystate, y, rhs, time); // damping buffer zone

  return;

} // add_RHSSourceTerms



void apply_physicalSource(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs, double time) {

  int num_samples = mystate->num_samples;

  // simulation
  if (mystate->simulation == "CASE_2DJET_WITH_A_HARMONIC_SOURCE")
    {
      // source parameters
      double amplitude = 0.001; // 0.001 kg/m/s^3
      double frequency = 76.0; // 76 rad/s
      double Bx = 0.04 * log(2.0); // 0.04 ln(2) m^{-2}
      double By = 0.32 * log(2.0); // 0.32 ln(2) m^{-2}

      // reference quantities
      double b = 1.3; // reference length, b = 1.3 m
      double R = 287.0; // gas constant, m^2/s^2/K
      double gamma_specificheat = mystate->gamma_specificheat;
      double rhoinfty = 1.184; // ambient density, 1.184 kg/m^3
      double Tinfty = 300.0; // ambient temperature, 300 K
      double cinfty = sqrt(gamma_specificheat * R * Tinfty); // ambient speed of sound

      // non-dimensionalize
      amplitude /= rhoinfty * pow(cinfty, 2) / (b / cinfty); // since it is a source to the pressure equation
      frequency *= math_constants::inverse_twopi * (b / cinfty);
      Bx *= pow(b, 2);
      By *= pow(b, 2);

      double angular_frequency = frequency * math_constants::twopi;
      double amp_cos_wt = amplitude * cos(angular_frequency * time);
      for (int l0 = 0; l0 < num_samples; l0++) {

        double xloc = mygrid->cell[l0].xyz[XDIR];
        double yloc = mygrid->cell[l0].xyz[YDIR];
        (rhs[IVAR_P])[l0] += exp(-(Bx * xloc*xloc + By * yloc*yloc)) * amp_cos_wt;

      } // l0
    }

  else if (mystate->simulation == "CASE_SCATTERING_TWOCYLINDER")
    {
      // source parameters
      // note that the source term defined in the category-2 NASA/CAA benchmark problem
      // (see its page 10 for equation 5) is for the energy equation as stated there
      // however, pressure equation is solved here and thus, a factor of gamma - 1 should 
      // be multiplied to be dimensionally consistent (since this is linear, it can be 
      // done at a post-processing stage of resulting pressure fluctuations though)
      double multi_factor = mystate->gamma_specificheat - 1.0;
      double exp_factor = -log(2.0) / (0.2 * 0.2);
      double sin_wt = sin( 4.0 * math_constants::twopi * time);
      double x0 = 0.0, y0 = 0.0;

      // as discussed by Manoha, Redonnet, Guenanff, and Terracol in the category-2 
      // NASA/CAA benchmark problem (see its page 277 for equation 6), the current 
      // simulation also suffers from large-amplitude wavefront generated at time 0 
      // and reflected back and forth within the computational domain
      // it turned out that this initial transient significantly corrupts simulation 
      // results; thus, the pressure source term is ramped up in time for its strength 
      // following their suggestions for 16 oscillation periods
      double ramp_factor = 1.0;
      double ramp_time = 16.0 / 4.0;
      if (time < ramp_time)
        ramp_factor = pow( sin(math_constants::pi / 2.0 * time / ramp_time), 2);

      for (int l0 = 0; l0 < num_samples; l0++) {

        double xs = mygrid->cell[l0].xyz[XDIR] - x0;
        double ys = mygrid->cell[l0].xyz[YDIR] - y0;
        (rhs[IVAR_P])[l0] += ramp_factor * multi_factor * exp(exp_factor * (xs*xs + ys*ys)) * sin_wt;

      } // l0
    }

  else if (mystate->simulation == "CASE_TANNA_TPN49")
    {
    }

  else if (mystate->simulation == "CASE_KBK_COMBUSTOR")
    {
    }

  else if (mystate->simulation == "CASE_LINEAR_NOZZLE")
    {
      // non-dimensionalize
      double amplitude = 1e-3;
      double frequency = 1.0 / 0.5042;
      double Bx = 20.0;
      double By = 20.0;

      double angular_frequency = frequency * math_constants::twopi;
      double amp_sin_wt = amplitude * sin(angular_frequency * time);
      for (int l0 = 0; l0 < num_samples; l0++) {
      
        double xloc = mygrid->cell[l0].xyz[XDIR] + 0.2;
        double rloc = mygrid->cell[l0].xyz[RDIR] - 0.0;
        //(rhs[IVAR_P])[l0] += exp(-(Bx * xloc*xloc + By * rloc*rloc)) * amp_sin_wt;
        (rhs[IVAR_P])[l0] += exp(-(Bx * xloc*xloc)) * amp_sin_wt;
      
      } // l0
    }

  else
    // do nothing since simulations may not have any source

  return;

} // apply_physicalSource



void apply_bufferZone(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs, double time) {

  int num_bufferZones = mygrid->num_bufferZones;

//  if ( time < 72.5) {
//    if ( mygrid->id_parent == 2 || mygrid->id_parent == 3 || mygrid->id_parent == 4 ) {
//      for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
//        for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
//          for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
//            int l0 = mygrid->idx1D(i, j, k);
//            if (mygrid->cell[l0].xyz[XDIR] > 10.0) {
//              double buffer_constant = 10.0;
//              buffer_constant = -buffer_constant/5.0 * (time - 67.5) + buffer_constant;
//              double buffer_strength = buffer_constant * pow((mygrid->cell[l0].xyz[XDIR] - 10.0)/90.0, 2);
//              for (int ivar = FIRST; ivar < mystate->num_vars_sol; ivar++) {
//                double relaxation = -buffer_strength * ((y[ivar])[l0] - mystate->sol_ref[ivar]); // relax to the prescribed ambient state
//                (rhs[ivar])[l0] += relaxation;
//              } // isol
//            } // mygrid->cell[l0].xyz[XDIR]
//          } // i
//        } // j
//      } // k
//    } // mygrid->id_parent
//  } // time

  for (int ibuffer = FIRST; ibuffer < num_bufferZones; ibuffer++) {

    Geometry::StructuredBufferZone *bufferZone_cur = &(mygrid->bufferZone[ibuffer]);

    switch ( bufferZone_cur->which_model ) {
    case SPONGE_FREUND_AMBIENT:

      apply_bufferZone_Freund_ambient(mygrid, mystate, y, rhs, time, bufferZone_cur);

      break;

    case SPONGE_FREUND_DIRICHLET:

      apply_bufferZone_Freund_Dirichlet(mygrid, mystate, y, rhs, time, bufferZone_cur);

      break;

    case SPONGE_FREUND_HARMONICWAVE:

      apply_bufferZone_Freund_harmonicWave(myinput, mygrid, mystate, y, rhs, time, bufferZone_cur);

      break;

    default:

      mpi::graceful_exit("Unknown type of buffer zone.");

      break;

    } // bufferZone_cur->which_model
  } // ibuffer

  return;

} // apply_bufferZone



void apply_bufferZone_Freund_ambient(Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs, double time, Geometry::StructuredBufferZone *mybufferZone) {

  for (int k = mybufferZone->is[ZETA]; k <= mybufferZone->ie[ZETA]; k++) {
    int k_in_grid = k - mybufferZone->is[ZETA] + mybufferZone->is_in_parent[ZETA];

    for (int j = mybufferZone->is[ETA]; j <= mybufferZone->ie[ETA]; j++) {
      int j_in_grid = j - mybufferZone->is[ETA] + mybufferZone->is_in_parent[ETA];

      for (int i = mybufferZone->is[XI]; i <= mybufferZone->ie[XI]; i++) {
        int i_in_grid = i - mybufferZone->is[XI] + mybufferZone->is_in_parent[XI];

        int lb = mybufferZone->idx1D(i, j, k);
        int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

        double buffer_strength = (mybufferZone->buffer_strength)[lb];
        for (int ivar = FIRST; ivar < mystate->num_vars_sol; ivar++) {

          double relaxation = -buffer_strength * ((y[ivar])[l0] - mystate->sol_ref[ivar]); // relax to the prescribed ambient state
          (rhs[ivar])[l0] += relaxation;

        } // isol
      } // i
    } // j
  } // k

  return;

} // apply_bufferZone_Freund_ambient



void apply_bufferZone_Freund_Dirichlet(Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs, double time, Geometry::StructuredBufferZone *mybufferZone) {

  // buffer zone is an object of a class derived from StructuredBoundaryCondition
  // thus, it can be also subject to a Dirichlet condition taking StructuredBoundaryCondition
  // as a first argument of mystate->prescribe_on_boundary_solution
  // the current physical model determines what kind of Dirichlet condition is applied

  int num_vars = mystate->num_vars_sol;
  int num_samples = mystate->num_samples;

  double **myboundarydata = new double *[num_vars];
  for (int ivar = 0; ivar < num_vars; ivar++)
    myboundarydata[ivar] = new double[mybufferZone->num_cells];

  for (int kb = mybufferZone->is[ZETA]; kb <= mybufferZone->ie[ZETA]; kb++)
    for (int jb = mybufferZone->is[ETA]; jb <= mybufferZone->ie[ETA]; jb++)
      for (int ib = mybufferZone->is[XI]; ib <= mybufferZone->ie[XI]; ib++) {

        int lb = mybufferZone->idx1D(ib, jb, kb);

        for (int ivar = 0; ivar < num_vars; ivar++)
          (myboundarydata[ivar])[lb] = DUMMY_DOUBLE; // some dummy value

      } // ib

  // the reference state is read from the Dirichlet boundary condition
  mystate->prescribe_on_boundary_solution(mybufferZone, mygrid, time, myboundarydata);

  for (int k = mybufferZone->is[ZETA]; k <= mybufferZone->ie[ZETA]; k++) {
    int k_in_grid = k - mybufferZone->is[ZETA] + mybufferZone->is_in_parent[ZETA];

    for (int j = mybufferZone->is[ETA]; j <= mybufferZone->ie[ETA]; j++) {
      int j_in_grid = j - mybufferZone->is[ETA] + mybufferZone->is_in_parent[ETA];

      for (int i = mybufferZone->is[XI]; i <= mybufferZone->ie[XI]; i++) {
        int i_in_grid = i - mybufferZone->is[XI] + mybufferZone->is_in_parent[XI];

        int lb = mybufferZone->idx1D(i, j, k);
        int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

        double buffer_strength = (mybufferZone->buffer_strength)[lb];
        for (int ivar = FIRST; ivar < num_vars; ivar++) {

          double relaxation = -buffer_strength * ((y[ivar])[l0] - (myboundarydata[ivar])[lb]);
          (rhs[ivar])[l0] += relaxation;

        } // isol
      } // i
    } // j
  } // k

  DEALLOCATE_2DPTR(myboundarydata, num_vars);

  return;

} // apply_bufferZone_Freund_Dirichlet



void apply_bufferZone_Freund_harmonicWave(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs, double time, Geometry::StructuredBufferZone *mybufferZone) {

  // buffer zone is an object of a class derived from StructuredBoundaryCondition
  // thus, it can be also subject to a Dirichlet condition taking StructuredBoundaryCondition
  // as a function argument

  int num_vars = mystate->num_vars_sol;
  int num_samples = mystate->num_samples;

  double **myboundarydata = new double *[num_vars];
  for (int ivar = 0; ivar < num_vars; ivar++)
    myboundarydata[ivar] = new double[mybufferZone->num_cells];

  for (int kb = mybufferZone->is[ZETA]; kb <= mybufferZone->ie[ZETA]; kb++)
    for (int jb = mybufferZone->is[ETA]; jb <= mybufferZone->ie[ETA]; jb++)
      for (int ib = mybufferZone->is[XI]; ib <= mybufferZone->ie[XI]; ib++) {

        int lb = mybufferZone->idx1D(ib, jb, kb);

        for (int ivar = 0; ivar < num_vars; ivar++)
          (myboundarydata[ivar])[lb] = DUMMY_DOUBLE; // some dummy value

      } // ib

  // the reference state is read from a time-harmonic wave solution
  // note that the time-harmonic wave parameters (e.g. amplitude, wavenumber, propagation direction) 
  // are specified in the input file
  bc::bc_dirichlet_harmonicwave(mybufferZone, mygrid, mystate, time, myboundarydata, myinput);

  for (int k = mybufferZone->is[ZETA]; k <= mybufferZone->ie[ZETA]; k++) {
    int k_in_grid = k - mybufferZone->is[ZETA] + mybufferZone->is_in_parent[ZETA];

    for (int j = mybufferZone->is[ETA]; j <= mybufferZone->ie[ETA]; j++) {
      int j_in_grid = j - mybufferZone->is[ETA] + mybufferZone->is_in_parent[ETA];

      for (int i = mybufferZone->is[XI]; i <= mybufferZone->ie[XI]; i++) {
        int i_in_grid = i - mybufferZone->is[XI] + mybufferZone->is_in_parent[XI];

        int lb = mybufferZone->idx1D(i, j, k);
        int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

        double buffer_strength = (mybufferZone->buffer_strength)[lb];
        for (int ivar = FIRST; ivar < num_vars; ivar++) {

          double relaxation = -buffer_strength * ((y[ivar])[l0] - (myboundarydata[ivar])[lb]);
          (rhs[ivar])[l0] += relaxation;

        } // isol
      } // i
    } // j
  } // k

  DEALLOCATE_2DPTR(myboundarydata, num_vars);

  return;

} // apply_bufferZone_Freund_harmonicWave

} // source
