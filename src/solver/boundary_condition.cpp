#include <cmath>

#include "boundary_condition.h"

namespace bc {

int num_vars;
int *varIndex_2update;

int num_of_boundaryCells;
int size_of_boundaryStencil;
double *stencil_boundary[MAX_ORDER_ACCURACY / 2];



void initialize_bc(int num_vars_in, int num_of_boundaryCells_in, int size_of_boundaryStencil_in, double *stencil_boundary_in[]) {

  num_vars = num_vars_in;

  // for some boundary conditions, it is sometimes desired not to update every solution variable
  // at the same time (so that different boundary condition routines can be re-used inside of the code)
  // this varIndex_2update makes that possible by storing the indices of variables to be updated
  // e.g. if you want to apply a boundary condition only to the 1st and 3rd variables, 
  //      varIndex_2update should be initialized by {0, 2} later at step 3
  //      also, make sure then to pass 2 instead of num_vars to a corresponding boundary 
  //      condition wrapper
  varIndex_2update = new int[num_vars];

  num_of_boundaryCells = num_of_boundaryCells_in;
  size_of_boundaryStencil = size_of_boundaryStencil_in;
  for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) {

    stencil_boundary[iboundary] = new double[size_of_boundaryStencil];

    for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++)
      (stencil_boundary[iboundary])[istencil] = (stencil_boundary_in[iboundary])[istencil];

  } // iboundary

  return;

} // initialize_bc



void apply_BC(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, int num_vars_in, double **y, double time) {

//  int irank_global_2lookat = 0;
//  int ivar_2lookat = IVAR_P;
//  if (mpi::irank == irank_global_2lookat) {
//
//    for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++)
//      for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++)
//        for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
//
//          int l0 = mygrid->idx1D(i, j, k);
//
//          (y[ivar_2lookat])[l0] = static_cast<double>(l0);
//          std::cout << "Rank: " << mpi::irank << ", ijk: " << i << " " << j << " " << k << ", val (before): " << (y[ivar_2lookat])[l0] << std::endl;
//
//        } // i
//
//  } // mpi::irank

  if (num_vars != num_vars_in)
    mpi::graceful_exit("The number of variables requiring boundary conditions becomes different than initialization.");

  int counter;

  for (int iboundary = 0; iboundary < mygrid->num_boundaryCondition_nonperiodic; iboundary++) {

    Geometry::StructuredBoundaryCondition *myboundary = &(mygrid->boundaryCondition[iboundary]);

    // step 1: ensure the boundary is subject to boundary condition
    //         boundaries doing overset-grid interpolation are bypassed
    if (myboundary->which_boundary != BOUNDARY_BC)
      continue;

    // step 2: temporary storage for boundary data
    double **data_boundary;
    ALLOCATE2D_DOUBLE(data_boundary, num_vars, myboundary->num_cells);
    for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++)
      for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++)
        for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {

          int lb = myboundary->idx1D(ib, jb, kb);

          for (int ivar = 0; ivar < num_vars; ivar++)
            (data_boundary[ivar])[lb] = -DUMMY_DOUBLE * static_cast<double>(myboundary->which_dir + 1); // some dummy value

        } // ib

    // step 3: provide boundary data
    for (int ivar = 0; ivar < num_vars; ivar++)
      varIndex_2update[ivar] = ivar; // by default, every variable is updated by a single type of boundary condition
    //
    switch( myboundary->which_model ) {
    case BC_DIRICHLET:

      mystate->prescribe_on_boundary_solution(myboundary, mygrid, time, data_boundary);

      break;

    case BC_DIRICHLET_ALLZERO:

      bc_dirichlet_allzero(myboundary, mygrid, time, data_boundary);

      break;

    case BC_DIRICHLET_HARMONICWAVE:

      bc_dirichlet_harmonicwave(myboundary, mygrid, time, data_boundary, myinput);

      break;

    case BC_NEUMANN:

      bc_neumann(myboundary, mygrid, time, data_boundary, num_vars, y);

      break;

    case BC_WALL_SLIP_KINEMATIC:

      // override varIndex_2update so that scalar variables are subject to the Neumann boundary condition
      for (int ivar = 0; ivar < num_vars; ivar++)
        varIndex_2update[ivar] = NONE;
      counter = 0;
      for (int ivar = 0; ivar < num_vars; ivar++)
        if (ivar != IVAR_UX && ivar != IVAR_UY && ivar != IVAR_UZ) // non-velocity variables are retained
          varIndex_2update[counter++] = ivar;

      // scalar variables are updated first by the Neumann boundary condition
      // note that the number of variables passed is decreased by DIM_MAX
      bc_neumann(myboundary, mygrid, time, data_boundary, num_vars - DIM_MAX, y);

      // update the remaining velocity components
      bc_wall_slip_kinematic(myboundary, mygrid, mystate, time, data_boundary, y);

      break;

    case BC_CENTERLINE_CART_NORM2X: // centerline in the Cartesian coordinates is perpendicular to a direction (either x, y, or z)
    case BC_CENTERLINE_CART_NORM2Y: // the velocity in that direction is zero on the centerline, while the other variables are subject to 
    case BC_CENTERLINE_CART_NORM2Z: // the Neumann boundary condition (symmetry)

      bc_neumann(myboundary, mygrid, time, data_boundary, num_vars, y);
      bc_centerline(myboundary, mygrid, time, data_boundary, CARTESIAN, myboundary->which_model - BC_CENTERLINE_CART_NORM2X);

      break;

    case BC_CENTERLINE_AXISYM: // the axisymmetric formulation assumes the cylindrical coorindates (x, r, \theta)
                               // thus, the value of the second direction is always zero at the axis of symmetry (i.e. r = 0)
                               // then, there is no need to branch into different cases like above Cartesian configurations

      bc_neumann(myboundary, mygrid, time, data_boundary, num_vars, y);
      bc_centerline(myboundary, mygrid, time, data_boundary, AXISYMMETRIC, RDIR);

      break;

    default:

      mpi::graceful_exit("Unknown type of boundary condition.");

      break;

    } // myboundary->which_model

    // step 4: update solution
    for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
      int k_in_grid = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

      for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
        int j_in_grid = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

        for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
          int i_in_grid = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

          int lb = myboundary->idx1D(ib, jb, kb);
          int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

          if (mygrid->cell[l0].iblank != BLANKED) // hole points are bypassed
            for (int ivar = 0; ivar < num_vars; ivar++)
              (y[ivar])[l0] = (data_boundary[ivar])[lb];

        } // ib
      } // jb
    } // kb

    // step 5: clean up
    DEALLOCATE_2DPTR(data_boundary, num_vars);

  } // iboundary

//  if (mpi::irank == irank_global_2lookat) {
//
//    for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++)
//      for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++)
//        for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
//
//          int l0 = mygrid->idx1D(i, j, k);
//
//          std::cout << "Rank: " << mpi::irank << ", ijk: " << i << " " << j << " " << k << ", val (after): " << (y[ivar_2lookat])[l0] << std::endl;
//
//        } // i
//
//  } // mpi::irank

  return;

} // apply_BC



void bc_dirichlet_allzero(Geometry::StructuredBoundaryCondition *myboundary, Geometry::StructuredGrid *mygrid, double time, double **myboundarydata) {

  for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++)
    for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++)
      for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {

        int lb = myboundary->idx1D(ib, jb, kb);

        for (int ivar = 0; ivar < num_vars; ivar++)
          (myboundarydata[ivar])[lb] = 0.0;

      } // ib

  return;

} // bc_dirichlet_allzero



void bc_dirichlet_harmonicwave(Geometry::StructuredBoundaryCondition *myboundary, Geometry::StructuredGrid *mygrid, double time, double **myboundarydata, UserInput *myinput) {

  int idx_in_grid[DIM_MAX];

  // ambient quantities
  double rho_0 = 1.0; // non-dimensional ambient density
  double p_0 = 1.0 / myinput->gamma_specificheat; // non-dimensional ambient pressure
  double c_0 = 1.0; // non-dimensional ambient speed of sound

  // time-harmonic wave parameters
  std::string waveType = myinput->harmonicWave_waveType;
  std::string waveForm = myinput->harmonicWave_waveForm;
  int idir_propagation = myinput->harmonicWave_idir_propagation;
  double amplitude = myinput->harmonicWave_amplitude;
  double wavelength = myinput->harmonicWave_wavelength;
  double wavenumber = math_constants::twopi / wavelength;

  double loc_propagation, loc_transverse[DIM_MAX-1];
  double pressure_fluctuation, velocity_fluctuation, entropy_fluctuation;

  for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
    idx_in_grid[ZETA] = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

    for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
      idx_in_grid[ETA] = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

      for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
        idx_in_grid[XI] = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

        int lb = myboundary->idx1D(ib, jb, kb);
        int l0 = mygrid->idx1D(idx_in_grid[XI], idx_in_grid[ETA], idx_in_grid[ZETA]);

        for (int ivar = 0; ivar < num_vars; ivar++)
          (myboundarydata[ivar])[lb] = 0.0;

        if (mygrid->cell[l0].iblank != BLANKED) { // hole points are bypassed

          loc_propagation = mygrid->cell[l0].xyz[idir_propagation];
          loc_transverse[FIRST] = mygrid->cell[l0].xyz[dir_other[idir_propagation][FIRST]];
          loc_transverse[SECOND] = mygrid->cell[l0].xyz[dir_other[idir_propagation][SECOND]];
          if (waveForm == "WAVEFORM_PLANE")
            if (waveType == "WAVE_ACOUSTIC") {

              pressure_fluctuation = amplitude * p_0 * sin(wavenumber * (loc_propagation - c_0 * time));
              velocity_fluctuation = pressure_fluctuation;
              entropy_fluctuation = 0.0;

            } // waveType
            else
              mpi::graceful_exit("Unknown type of wave for this wave form.");

          else if (waveForm == "WAVEFORM_HOMOGENEOUS")
            if (waveType == "WAVE_ACOUSTIC") {

              pressure_fluctuation = amplitude * p_0 * sin(wavenumber * c_0 * time);
              velocity_fluctuation = pressure_fluctuation;
              entropy_fluctuation = 0.0;

            } // waveType
            else if (waveType == "WAVE_ENTROPY") {

              pressure_fluctuation = 0.0;
              velocity_fluctuation = 0.0;
              entropy_fluctuation = amplitude * (1.0 + sin(wavenumber * c_0 * time)); // only non-negative entropy fluctuations are considered

            } // waveType
            else
              mpi::graceful_exit("Unknown type of wave for this wave form.");

          else if (waveForm == "WAVEFORM_GAUSSIAN_ROUNDJET")
            if (waveType == "WAVE_ENTROPY") {

              pressure_fluctuation = 0.0;
              velocity_fluctuation = 0.0;
              entropy_fluctuation = amplitude * (1.0 + sin(wavenumber * c_0 * time)); // only non-negative entropy fluctuations are considered
              //
              double fac_Gaussian = -log(2.0) / pow(myinput->harmonicWave_halfWidth, 2);
              entropy_fluctuation *= exp(fac_Gaussian * (pow(loc_transverse[FIRST], 2) + 
                                                         pow(loc_transverse[SECOND], 2)));

            } // waveType
            else
              mpi::graceful_exit("Unknown type of wave for this wave form.");

          else
            mpi::graceful_exit("Unknown type of wave form.");

          if (myinput->model_pde == "LINEAR_ACOUSTICS") {

            (myboundarydata[IVAR_RHO])[lb] = pressure_fluctuation / pow(c_0, 2.0);
            (myboundarydata[IVAR_UX + idir_propagation])[lb] = velocity_fluctuation / (rho_0 * c_0);
            (myboundarydata[IVAR_P])[lb] = pressure_fluctuation;

          } // myinput->model_pde
          else if (myinput->model_pde == "LEE" ||
                   myinput->model_pde == "LEE_SCALAR") {

            (myboundarydata[IVAR_S])[lb] = entropy_fluctuation;
            (myboundarydata[IVAR_UX + idir_propagation])[lb] = velocity_fluctuation / (rho_0 * c_0);
            (myboundarydata[IVAR_P])[lb] = pressure_fluctuation;

          } // myinput->model_pde
          else
            mpi::graceful_exit("The Dirichlet boundary enforcing a time-harmonic wave is not implemented for this physical model.");

        } // mygrid->cell[l0].iblank
      } // ib
    } // jb
  } // kb

  return;

} // bc_dirichlet_harmonicwave



void bc_neumann(Geometry::StructuredBoundaryCondition *myboundary, Geometry::StructuredGrid *mygrid, double time, double **myboundarydata, int num_vars_in, double **y) {

  int idir_boundary = myboundary->which_dir; // XI or ETA or ZETA
  int iend_boundary = myboundary->which_end; // LEFT or RIGHT
  int idx_in_grid[DIM_MAX];
  int multi_fac;

  switch ( idir_boundary ) {
  case XI:
    multi_fac = 1;

    break;

  case ETA:
    multi_fac = mygrid->nXi;

    break;

  case ZETA:
    multi_fac = mygrid->nXi * mygrid->nEta;

    break;

  } // idir_boundary

  if (iend_boundary == LEFT) {
      for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
        idx_in_grid[ZETA] = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

        for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
          idx_in_grid[ETA] = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

          for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
            idx_in_grid[XI] = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

            int lb = myboundary->idx1D(ib, jb, kb);
            int l0 = mygrid->idx1D(idx_in_grid[XI], idx_in_grid[ETA], idx_in_grid[ZETA]);

            if (mygrid->cell[l0].iblank != BLANKED) { // hole points are bypassed

              for (int ivar = 0; ivar < num_vars_in; ivar++) {
                int idx_in_var = varIndex_2update[ivar];

                (myboundarydata[idx_in_var])[lb] = 0.0;
                for (int istencil = SECOND; istencil < size_of_boundaryStencil; istencil++) {

                  int l1 = l0 + multi_fac * istencil;
                  (myboundarydata[idx_in_var])[lb] -= stencil_boundary[FIRST][istencil] * (y[idx_in_var])[l1];

                } // istencil
                (myboundarydata[idx_in_var])[lb] /= stencil_boundary[FIRST][FIRST];

              } // ivar
            } // mygrid->cell[l0].iblank
          } // ib
        } // jb
      } // kb

  } // iend_boundary
  else if (iend_boundary == RIGHT) {
    for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
      idx_in_grid[ZETA] = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

      for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
        idx_in_grid[ETA] = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

        for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
          idx_in_grid[XI] = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

          int lb = myboundary->idx1D(ib, jb, kb);
          int l0 = mygrid->idx1D(idx_in_grid[XI], idx_in_grid[ETA], idx_in_grid[ZETA]);

          if (mygrid->cell[l0].iblank != BLANKED) { // hole points are bypassed

            for (int ivar = 0; ivar < num_vars_in; ivar++) {
              int idx_in_var = varIndex_2update[ivar];

              (myboundarydata[idx_in_var])[lb] = 0.0;
              for (int istencil = SECOND; istencil < size_of_boundaryStencil; istencil++) {

                int l1 = l0 - multi_fac * istencil;
                (myboundarydata[idx_in_var])[lb] += stencil_boundary[FIRST][istencil] * (y[idx_in_var])[l1];

              } // istencil
              (myboundarydata[idx_in_var])[lb] /= -stencil_boundary[FIRST][FIRST];

            } // ivar
          } // mygrid->cell[l0].iblank
        } // ib
      } // jb
    } // kb

  } // iend_boundary

  return;

} // bc_neumann



void bc_wall_slip_kinematic(Geometry::StructuredBoundaryCondition *myboundary, Geometry::StructuredGrid *mygrid, State *mystate, double time, double **myboundarydata, double **y) {

  int idir_boundary = myboundary->which_dir; // XI or ETA or ZETA

  double velocity_Cartesian[DIM_MAX];
  double velocity_contravariant[DIM_MAX];

  for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
    int k_in_grid = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

    for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
      int j_in_grid = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

      for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
        int i_in_grid = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

        int lb = myboundary->idx1D(ib, jb, kb);
        int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

        for (int idir = XDIR; idir < DIM_MAX; idir++)
          velocity_Cartesian[idir] = 0.0;

        if (mygrid->cell[l0].iblank != BLANKED) { // hole points are bypassed

          // compute contravariant velocity components
          for (int idir = XDIR; idir < DIM_MAX; idir++)
            velocity_Cartesian[idir] = (y[IVAR_UX + idir])[l0];
          mystate->to_contravariant_velocity(mygrid->cell[l0].metrics, velocity_Cartesian, velocity_contravariant);

          // impose the non-penetration condition
          velocity_contravariant[idir_boundary] = 0.0;

          // back to Cartesian velocity components
          mystate->to_Cartesian_velocity(mygrid->cell[l0].metricsInverse, mygrid->cell[l0].Jac, velocity_contravariant, velocity_Cartesian);

        } // mygrid->cell[l0].iblank

        for (int idir = XDIR; idir < DIM_MAX; idir++)
          (myboundarydata[IVAR_UX + idir])[lb] = velocity_Cartesian[idir];

      } // ib
    } // jb
  } // kb

  return;

} // bc_wall_slip_kinematic



void bc_centerline(Geometry::StructuredBoundaryCondition *myboundary, Geometry::StructuredGrid *mygrid, double time, double **myboundarydata, int which_coord_system, int which_dir) {

  int which_vel;

  // all this subroutine does is to locate a velocity component perpendicular to the centerline and to zero it out
  // this is only true for the Cartesian coordinates where the normal velocity component is anti-symmetric
  // for the cylindrical coordinates, all the velocity components other than the one in the centerline direction should be zero

  switch ( which_coord_system ) {
  case CARTESIAN:

    which_vel = IVAR_UX + which_dir;

    for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
      int k_in_grid = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

      for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
        int j_in_grid = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

        for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
          int i_in_grid = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

          int lb = myboundary->idx1D(ib, jb, kb);
          int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

          if (mygrid->cell[l0].iblank != BLANKED) { // hole points are bypassed

            (myboundarydata[which_vel])[lb] = 0.0;

          } // mygrid->cell[l0].iblank
        } // ib
      } // jb
    } // kb

    break;

  case AXISYMMETRIC: // since axisymmetry requires 2-D, no need to attempt to zero out 
                     // the velocity components along the other two directions than 
                     // the axial direction (i.e. x); it is sufficient to work with 
                     // the radial component

    which_vel = IVAR_UR;

    for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
      int k_in_grid = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

      for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
        int j_in_grid = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

        for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
          int i_in_grid = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

          int lb = myboundary->idx1D(ib, jb, kb);
          int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

          if (mygrid->cell[l0].iblank != BLANKED) { // hole points are bypassed

            (myboundarydata[which_vel])[lb] = 0.0;

          } // mygrid->cell[l0].iblank
        } // ib
      } // jb
    } // kb

    break;

  default:

    mpi::graceful_exit("Unknown coordinate system to apply a centerline boundary condition.");

    break;

  } // which_coord_system

  return;

} // bc_centerline

} // bc
