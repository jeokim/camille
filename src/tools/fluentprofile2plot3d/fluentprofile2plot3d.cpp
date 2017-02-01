#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "../../core/simulation.h"
#include "../../io/fluent.h"
#include "../../solver/spatial_discretization.h"

int main(int argc, char * argv[]) {

  // global variables and arrays
  int num_dim_source;
  int num_vars_source;
  int num_cvs_FLUENT = 0;
  double **xyz_source;
  double **var_source;
  std::string str_stdout;
  std::stringstream str_counter;

  // simplicity
  UserInput *myinput = &(simulation::myinput);
  Geometry::StructuredBlock *&block = simulation::block;
  Geometry::StructuredGrid *&mygrid = simulation::grid; // "Geometry::StructuredGrid *mygrid = simulation::grid;" doesn't work 
                                                        // since mygrid gets only a copy of "simulation::grid" although the syntax is 
                                                        // consistent (a pointer copied to a pointer); to access the original 
                                                        // "simulation::grid", it is necessary to first make the receiver "mygrid" as 
                                                        // a reference (by putting & in front of "mygrid"); then, to be consistent (a 
                                                        // pointer on LHS and a pointer on RHS), place * in its front
  State *mystate = &(simulation::state);



  // since the interpolation is done in parallel and filtering is applied, the 
  // entire initialization routine is called
  simulation::initialize(argc, argv);



  // some temporary hacks
  if (mpi::irank == 0)
    std::cout << std::endl;
  MESSAGE_STDOUT("The FLUENT--PLOT3D interpolation is currently limited to 2D-2D cases (axisymmetric case included).");
  MESSAGE_STDOUT("The FLUENT profile should contain absolute pressure, density, axial (or streamwise) velocity, radial (or transverse) velocity.");
  if (myinput->num_dim != 2)
    mpi::graceful_exit("Only 2D-2D interpolation is supported.");



  // read a FLUENT profile
  num_dim_source = myinput->num_dim_interpSource;
  num_vars_source = myinput->num_vars_interpSource;
  fluent::read_profile(myinput->file_profile_FLUENT, num_dim_source, num_vars_source, num_cvs_FLUENT, xyz_source, var_source);
  mpi::wait_allothers("A interpolation-source file is read and stored.");



  // re-scale the xyz locations of the source data
  double inv_scale_xyz = 1.0 / myinput->scale_xyz;
  for (int icv = 0; icv < num_cvs_FLUENT; icv++)
    for (int idim = XDIR; idim < num_dim_source; idim++)
      xyz_source[icv][idim] *= inv_scale_xyz;
  mpi::wait_allothers("The spatial locations of the interpolation-source data are re-scaled.");



  // post-process the solution to be interpolated
  double gamma;
  if (myinput->model_pde == "LEE") {

    // re-scale the solution to be interpolated
    MESSAGE_STDOUT("P, RHO, U_i are assumed to be in the interpolation-source file in that order.");
    double inv_scale_p = 1.0 / myinput->scale_p;
    double inv_scale_rho = 1.0 / myinput->scale_rho;
    double inv_scale_u = 1.0 / myinput->scale_u;
    for (int icv = 0; icv < num_cvs_FLUENT; icv++) {

      var_source[icv][FIRST] *= inv_scale_p; // P
      var_source[icv][SECOND] *= inv_scale_rho; // RHO
      var_source[icv][THIRD] *= inv_scale_u; // U_X
      var_source[icv][FOURTH] *= inv_scale_u; // U_Y

    } // icv
    mpi::wait_allothers("The solution to be interpolated is non-dimensionalized.");

    // convert {P, RHO, U_i} into {S, U_i, P}
    gamma = myinput->gamma_specificheat;
    for (int icv = 0; icv < num_cvs_FLUENT; icv++) {

      double entropy = (log(gamma * var_source[icv][FIRST]) - gamma * log(var_source[icv][SECOND])) / gamma;
      double pressure = var_source[icv][FIRST];
      double ux = var_source[icv][THIRD];
      double uy = var_source[icv][FOURTH];

      var_source[icv][FIRST] = entropy;
      var_source[icv][SECOND] = pressure;
      var_source[icv][THIRD] = ux;
      var_source[icv][FOURTH] = uy;

    } // icv
    mpi::wait_allothers("The solution to be interpolated is converted to be consistent with the current solution.");

  } // myinput->model_pde
  else
    mpi::graceful_exit("Unsupported physical model for the interpolation.");



  // do a nearest-neighbor interpolation
  double xyz[DIM_MAX];
  int icv_distance_min;
  double distance, distance_min;
  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int l0 = mygrid->idx1D(i, j, k);

        // the xyz locations of a current point (e.g. PLOT3D) to be interpolated
        for (int idim = XDIR; idim < DIM_MAX; idim++)
          xyz[idim] = mygrid->cell[l0].xyz[idim];

        // find the minimum distance to the FLUENT points by a brute force
        icv_distance_min = NONE;
        distance_min = DUMMY_LARGE;
        for (int icv = 0; icv < num_cvs_FLUENT; icv++) {

          distance = 0.0;
          for (int idim = XDIR; idim < num_dim_source; idim++)
            distance += pow(xyz[idim] - xyz_source[icv][idim], 2); // sqrt is not applied since only the closest-point's index is of interest

          if (distance < distance_min) {

            icv_distance_min = icv;
            distance_min = distance;

          } // distance
        } // icv
        assert( icv_distance_min != NONE );
        assert( distance_min != DUMMY_LARGE );
        distance_min = sqrt( distance_min );

        // get the solution of the nearest neighbor
        (mystate->sol[IVAR_S])[l0] = var_source[icv_distance_min][FIRST];
        (mystate->sol[IVAR_UX])[l0] = var_source[icv_distance_min][THIRD];
        (mystate->sol[IVAR_UY])[l0] = var_source[icv_distance_min][FOURTH];
        (mystate->sol[IVAR_P])[l0] = var_source[icv_distance_min][SECOND];

      } // i
    } // j
  } // k
  mpi::wait_allothers("The nearest neighbor interpolation is completed.");



  // apply a filter
  for (int ifilter = 0; ifilter < myinput->num_filters_interpolation; ifilter++) {

    spatial::apply_filter(mystate->num_vars_sol, mystate->num_samples, mystate->sol, mygrid);
    if (mpi::irank == 0)
      if (ifilter == 0)
        std::cout << ">>> The filter is applied " << std::setw(6) << ifilter + 1 << " time." << std::endl;
      else
        std::cout << ">>> The filter is applied " << std::setw(6) << ifilter + 1 << " times." << std::endl;

  } // ifilter



  // write the result
  mystate->time_step = 0;
  mystate->time_sol = 0.0;
  if (myinput->interpolate_into_which_PLOT3D == "PLOT3D_SOLUTION") {

    io::write_solution(myinput, mygrid, block, mystate);
    mpi::wait_allothers("The interpolated solution is written; check the regular solution file written at the counter 0.");

  } // myinput->interpolate_into_which_PLOT3D
  else if (myinput->interpolate_into_which_PLOT3D == "PLOT3D_FUNCTION") {

    for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
      for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
        for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
          int l0 = mygrid->idx1D(i, j, k);

          for (int ivar = FIRST; ivar < mystate->num_vars_sol; ivar++)
            (mystate->sol_mean[ivar])[l0] = (mystate->sol[ivar])[l0];

        } // i
      } // j
    } // k

    io::write_solution_mean(myinput, mygrid, block, mystate);
    mpi::wait_allothers("The interpolated mean state is written; check the solution file starting with 'mean_'.");

  } // myinput->interpolate_into_which_PLOT3D

  // clean up
  DEALLOCATE_2DPTR(xyz_source, num_cvs_FLUENT);
  DEALLOCATE_2DPTR(var_source, num_cvs_FLUENT);

  return 0;

} // main
