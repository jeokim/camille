#include <cmath>
#include <iostream>
#include <sstream>

#include "state.h"

void State::initialize_state(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // dimensionality
  this->num_dim = myinput->num_dim;

  // physical model
  this->model_pde = myinput->model_pde;

  // simulation
  this->simulation = myinput->simulation;

  // solution size
  this->num_samples = mygrid->num_ocells;
  this->num_cells_ghost = myinput->num_cells_ghost;

  // time information
  this->time_step = 0;
  this->time_step_lastrun = 0;
  this->time_sol = 0.0;
  this->time_sol_lastrun = 0.0;

  // set a number of variables
  this->num_vars_sol = myinput->num_vars_sol;
  this->num_vars_mean = myinput->num_vars_mean;
  this->num_vars_meanGradient = myinput->num_vars_meanGradient;
  this->num_vars_aux = myinput->num_vars_aux;

  // allocate
  this->sol = new double *[this->num_vars_sol];
  this->sol_old = new double *[this->num_vars_sol];
  for (int ivar = 0; ivar < this->num_vars_sol; ivar++) {

    this->sol[ivar] = new double[this->num_samples];
    this->sol_old[ivar] = new double[this->num_samples];

  } // ivar
  this->sol_ref = new double[this->num_vars_sol];
  if (this->num_vars_mean > 0) { // if base or mean states are used for the solver

    this->sol_mean = new double *[this->num_vars_mean];
    for (int ivar = 0; ivar < this->num_vars_mean; ivar++)
      this->sol_mean[ivar] = new double[this->num_samples];

  } // this->num_vars_mean
  if (this->num_vars_meanGradient > 0) { // if gradients of base or mean states are used for the solver

    this->sol_meanGradient = new double *[this->num_vars_mean * this->num_dim]; // it is a gradient vector; thus, multiply by a current dimension
    for (int ivar = 0; ivar < this->num_vars_mean * this->num_dim; ivar++)
      this->sol_meanGradient[ivar] = new double[this->num_samples];

  } // this->num_vars_meanGradient
  if (this->num_vars_aux > 0) { // if there are any auxiliary variables used for the solver

    this->sol_aux = new double *[this->num_vars_aux];
    for (int ivar = 0; ivar < this->num_vars_aux; ivar++)
      this->sol_aux[ivar] = new double[this->num_samples];

  } // this->num_vars_aux

  // initialize (generic)
  for (int ivar = 0; ivar < this->num_vars_sol; ivar++) {
    for (int l0 = 0; l0 < this->num_samples; l0++) {

      (this->sol[ivar])[l0] = 0.0;
      (this->sol_old[ivar])[l0] = 0.0;

    } // l0

    this->sol_ref[ivar] = 0.0;

  } // ivar
  if (this->num_vars_mean > 0) // if base or mean states are used for the solver
    for (int ivar = 0; ivar < this->num_vars_mean; ivar++)
      for (int l0 = 0; l0 < this->num_samples; l0++)
        (this->sol_mean[ivar])[l0] = 0.0;
  if (this->num_vars_meanGradient > 0) // if gradients of base or mean states are used for the solver
    for (int ivar = 0; ivar < this->num_vars_mean * this->num_dim; ivar++)
      for (int l0 = 0; l0 < this->num_samples; l0++)
        (this->sol_meanGradient[ivar])[l0] = 0.0;
  if (this->num_vars_aux > 0) // if there are any auxiliary variables used for the solver
    for (int ivar = 0; ivar < this->num_vars_aux; ivar++)
      for (int l0 = 0; l0 < this->num_samples; l0++)
        (this->sol_aux[ivar])[l0] = 0.0;

  this->gamma_specificheat = myinput->gamma_specificheat;

  // physical model
  if (this->model_pde == "LINEAR_ACOUSTICS")
    this->initialize_state_acoustics(myinput, mygrid);

  else if (this->model_pde == "LEE")
    this->initialize_state_linearizedEuler(myinput, mygrid);

  else if (this->model_pde == "LEE_SCALAR") {
    this->initialize_state_linearizedEuler(myinput, mygrid);
    this->initialize_state_linearizedEuler_scalar(myinput, mygrid);
  } // this->model_pde

  else if (this->model_pde == "LEE_MIXFRAC_CONSTGAMMA") {
    this->initialize_state_linearizedEuler(myinput, mygrid);
    this->initialize_state_linearizedEuler_scalar(myinput, mygrid);
    this->initialize_state_linearizedEuler_aux_composition(myinput, mygrid);
  } // this->model_pde

  else
    mpi::graceful_exit("PHYSICAL_MODEL " + this->model_pde + " is not implemented and cannot be initialized.");

  // simulation
  if (this->simulation == "CASE_PLANE_WAVE")
    this->initialize_state_acoustics_harmonicWave(myinput, mygrid);

  else if (this->simulation == "CASE_GAUSSIAN_PULSE")
    this->initialize_state_acoustics_GaussianPulse(myinput, mygrid);

  else if (this->simulation == "CASE_2DJET_WITH_A_HARMONIC_SOURCE")
    this->initialize_state_linearizedEuler_2Djet_with_a_harmonic_source(myinput, mygrid);

  else if (this->simulation == "CASE_SCATTERING_TWOCYLINDER")
    this->initialize_state_linearizedEuler_cylinderScattering(myinput, mygrid);

  else if (this->simulation == "CASE_TANNA_TPN49")
    this->initialize_state_linearizedEuler_TannaTPN49(myinput, mygrid);

  else if (this->simulation == "CASE_KBK_COMBUSTOR")
    this->initialize_state_linearizedEuler_KBKCombustor(myinput, mygrid);

  else if (this->simulation == "CASE_LINEAR_NOZZLE")
    this->initialize_state_linearizedEuler_linearNozzle(myinput, mygrid);

  else
    mpi::graceful_exit("SIMULATION " + this->simulation + " is not implemented and cannot be initialized.");

  this->compute_dependent_variables(this->sol);

  return;

} // State::initialize_state



void State::initialize_state_acoustics(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // linear acoustics
  //   solution variables: rho', velocity', p'
  //   auxiliary variables: none
  // see core/param.h to see the variable ordering

  // set names of variables
  this->name_vars = new std::string[this->num_vars_sol];
  this->name_vars[IVAR_RHO] = "RHO'";
  this->name_vars[IVAR_UX] = "UX'";
  this->name_vars[IVAR_UY] = "UY'";
  this->name_vars[IVAR_UZ] = "UZ'";
  this->name_vars[IVAR_P] = "P'";
  //
  if (myinput->axisym == TRUE) {

    this->name_vars[IVAR_UR] = "UR'";
    this->name_vars[IVAR_UTHETA] = "UTHETA'";

  } // myinput->axisym

  // ambient reference state
  for (int ivar = IVAR_RHO; ivar <= IVAR_P; ivar++)
    this->sol_ref[ivar] = 0.0; // ambient state of fluctuating variables is all zero

  // initialize
  for (int ivar = IVAR_RHO; ivar <= IVAR_P; ivar++)
    for (int l0 = 0; l0 < this->num_samples; l0++)
      (this->sol[ivar])[l0] = 0.0;

  return;

} // State::initialize_state_acoustics



void State::initialize_state_acoustics_harmonicWave(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

//  int idir_wave_propagation = XDIR;
  double time = 0.0;

  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
      int k_in_block = k - mygrid->iso[ZETA] + mygrid->iso_in_parent[ZETA];

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
        int j_in_block = j - mygrid->iso[ETA] + mygrid->iso_in_parent[ETA];

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
          int i_in_block = i - mygrid->iso[XI] + mygrid->iso_in_parent[XI];

        int l0 = mygrid->idx1D(i, j, k);

//        double *plane_wave = this->acoustics_plane_wave(time, mygrid->cell[l0].xyz[idir_wave_propagation], idir_wave_propagation);

//        for (int ivar = IVAR_RHO; ivar <= IVAR_P; ivar++)
//          (this->sol[ivar])[l0] = plane_wave[ivar];
//          (this->sol[ivar])[l0] = pow( mygrid->cell[l0].xyz[YDIR], 2*ivar);
//          (this->sol[ivar])[l0] = static_cast<double>(mpi::irank+1);

//        DEALLOCATE_1DPTR(plane_wave);

      } // i
    } // j
  } // k

  return;

} // State::initialize_state_acoustics_harmonicWave



double *State::acoustics_plane_wave(double time, double location, int idir_propagation) {

  // ambient quantities
  double rho_0 = 1.0;
  double p_0 = 1.0 / this->gamma_specificheat;
  double c_0 = 1.0; // ambient speed of sound
  double g_0 = 0.005; // relative magnitude of pressure fluctuation

  double wavelength = 4.0;
  double wavenumber = 2.0 * math_constants::pi / wavelength;

  double pressure_fluctuation = g_0 * p_0 * sin(wavenumber * (location - c_0 * time));

  double *sol_tmp = new double[this->num_vars_sol];
  sol_tmp[IVAR_RHO] = pressure_fluctuation / pow(c_0, 2.0);
  for (int ivel = IVAR_UX; ivel <= IVAR_UZ; ivel++)
    sol_tmp[ivel] = 0.0;
  sol_tmp[IVAR_UX + idir_propagation] = pressure_fluctuation / (rho_0 * c_0);
  sol_tmp[IVAR_P] = pressure_fluctuation;

  return sol_tmp;

} // State::acoustics_plane_wave



void State::initialize_state_acoustics_GaussianPulse(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  double rho_0 = 1.0; // ambient density
  double p_0 = 1.0 / this->gamma_specificheat; // ambient pressure (1/gamma)
  double c_0 = 1.0; // ambient speed of sound
  double g_0 = 0.005; // relative magnitude of pressure fluctuation
  //
  double x_0 = 15.0;
  double y_0 = 0.0;
  double z_0 = 0.0;
  double sigma = pow(3.0, 2);

  for (int k = mygrid->iso[ZETA]; k < mygrid->ieo[ZETA] + 1; k++) {
    for (int j = mygrid->iso[ETA]; j < mygrid->ieo[ETA] + 1; j++) {
      for (int i = mygrid->iso[XI]; i < mygrid->ieo[XI] + 1; i++) {

        int l0 = mygrid->idx1D(i, j, k);
        double radius = pow(mygrid->cell[l0].xyz[XDIR] - x_0, 2) + 
                        pow(mygrid->cell[l0].xyz[YDIR] - y_0, 2) + 
                        pow(mygrid->cell[l0].xyz[ZDIR] - z_0, 2);
//        double radius = pow(mygrid->cell[l0].xyz[XDIR] - x_0, 2);
        double pressure_fluctuation = g_0 * p_0 * exp(-sigma * radius);

        (this->sol[IVAR_RHO])[l0] = pressure_fluctuation / pow(c_0, 2.0);
        (this->sol[IVAR_UX])[l0] = 0.0;
        (this->sol[IVAR_UY])[l0] = 0.0;
        (this->sol[IVAR_UZ])[l0] = 0.0;
        (this->sol[IVAR_P])[l0] = pressure_fluctuation;

      } // i
    } // j
  } // k

  return;

} // State::initialize_state_acoustics_GaussianPulse



void State::initialize_state_linearizedEuler(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // linearized Euler model
  //   solution variables: s', velocity', p'
  //   auxiliary variables: rho', mean of rho, T', mean of T
  // see core/param.h to see the variable ordering

  // set names of variables
  this->name_vars = new std::string[this->num_vars_sol];
  this->name_vars[IVAR_S] = "S'";
  this->name_vars[IVAR_UX] = "UX'";
  this->name_vars[IVAR_UY] = "UY'";
  this->name_vars[IVAR_UZ] = "UZ'";
  this->name_vars[IVAR_P] = "P'";
  //
  this->name_vars_mean = new std::string[this->num_vars_mean];
  this->name_vars_mean[IVAR_S] = "Sbar";
  this->name_vars_mean[IVAR_UX] = "UXbar";
  this->name_vars_mean[IVAR_UY] = "UYbar";
  this->name_vars_mean[IVAR_UZ] = "UZbar";
  this->name_vars_mean[IVAR_P] = "Pbar";
  //
  this->name_vars_aux = new std::string[this->num_vars_aux];
  this->name_vars_aux[IAUX_RHO] = "RHO'";
  this->name_vars_aux[IAUX_RHO_MEAN] = "RHObar";
  this->name_vars_aux[IAUX_T] = "T'";
  this->name_vars_aux[IAUX_T_MEAN] = "Tbar";
  //
  if (myinput->axisym == TRUE) {

    this->name_vars[IVAR_UR] = "UR'";
    this->name_vars[IVAR_UTHETA] = "UTHETA'";
    //
    this->name_vars_mean[IVAR_UR] = "URbar";
    this->name_vars_mean[IVAR_UTHETA] = "UTHETAbar";

  } // myinput->axisym

  // ambient reference state
  for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
    this->sol_ref[ivar] = 0.0; // ambient state of fluctuating variables is all zero

  // initialize
  for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
    for (int l0 = 0; l0 < this->num_samples; l0++)
      (this->sol[ivar])[l0] = 0.0;
  for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
    for (int l0 = 0; l0 < this->num_samples; l0++)
      (this->sol_mean[ivar])[l0] = 0.0;
  for (int ivar = 0; ivar < this->num_vars_aux; ivar++)
    for (int l0 = 0; l0 < this->num_samples; l0++)
      (this->sol_aux[ivar])[l0] = 0.0;

  return;

} // State::initialize_state_linearizedEuler



void State::initialize_state_linearizedEuler_scalar(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // linearized Euler model
  //   solution variables: s', velocity', p'
  //   auxiliary variables: rho', mean of rho, T', mean of T
  // see core/param.h to see the variable ordering

  // only passive scalar variables are initialized

  // so that scalars are found right behind p'
  int ivar_shift = IVAR_P + 1;
  std::string varname = "Z";
  std::stringstream str_counter;

  // set names of variables
  for (int ivar = 0; ivar < myinput->num_scalar; ivar++) {
    str_counter << ivar; // scalar names go like Z0, Z1, Z2, ...

    // for more than one scalar, name them as Z0', Z1', ...
    // for a single scalar, use Z'
    if (myinput->num_scalar > 1) {
      this->name_vars[ivar_shift+ivar] = varname + str_counter.str() + "'";
      this->name_vars_mean[ivar_shift+ivar] = varname + str_counter.str() + "bar";
    } // myinput->num_scalar
    else {
      this->name_vars[ivar_shift+ivar] = varname + "'";
      this->name_vars_mean[ivar_shift+ivar] = varname + "bar";
    }
  } // ivar

  // ambient reference state
  for (int ivar = 0; ivar < myinput->num_scalar; ivar++)
    this->sol_ref[ivar_shift+ivar] = 0.0; // ambient state of fluctuating variables is all zero

  // initialize
  for (int ivar = 0; ivar < myinput->num_scalar; ivar++)
    for (int l0 = 0; l0 < this->num_samples; l0++)
      (this->sol[ivar_shift+ivar])[l0] = 0.0;
  for (int ivar = 0; ivar < myinput->num_scalar; ivar++)
    for (int l0 = 0; l0 < this->num_samples; l0++)
      (this->sol_mean[ivar_shift+ivar])[l0] = 0.0;

  return;

} // State::initialize_state_linearizedEuler_scalar



void State::initialize_state_linearizedEuler_aux_composition(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // linearized Euler model
  //   solution variables: s', velocity', p'
  //   auxiliary variables: rho', mean of rho, T', mean of T
  // see core/param.h to see the variable ordering

  // only additional auxiliary variables (c_p, dc_p/dZ, & \Psi) are updated

  // set names of variables
  this->name_vars_aux[IAUX_CP] = "CP";
  this->name_vars_aux[IAUX_DCPDZ] = "DCPDZ";
  this->name_vars_aux[IAUX_PSI] = "PSI";

  return;

} // State::initialize_state_linearizedEuler_aux



void State::initialize_state_linearizedEuler_2Djet_with_a_harmonic_source(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // harmonic acoustic source within a 2-D parallel jet
  // see the fourth computational aeroacoustics workshop on benchmark problems (2004)
  // also see Agarwal et al. (AIAA J. Vol. 42, No. 1, 2004)

  // reference length scale: b in eq. (18) (equivalent to the jet half width)
  // reference density scale: \rho_\infty
  // reference velocity scale: c_\infty
  // reference pressure scale: \rho_\infty c_\infty^2
  // reference temperature scale: (\gamma - 1) T_\infty
  // reference entropy scale: c_\infty^2 / ((\gamma - 1) T_\infty) = Cp

  // ambient state
  double R = 287.0; // gas constant, m^2/s^2/K
  double Cv = R / (this->gamma_specificheat - 1.0); // specific heat at a constant volume
  double Cp = Cv + R; // specific heat at a constant pressure, Cp - Cv = R
  double rhoinfty = 1.184; // ambient density, 1.184 kg/m^3
  double Tinfty = 300.0; // ambient temperature, 300 K
  double cinfty = sqrt(this->gamma_specificheat * R * Tinfty); // ambient speed of sound
  if (mpi::irank == 0) {

    std::cout << std::endl;
    std::cout << "Variables at the ambient state: " << std::endl;
    std::cout << "  Speed of sound: " << cinfty << " [m/s]" << std::endl;

  } // mpi::irank

  // base state
//  double pbar = 103330.0; // mean pressure, 103330 kg/m/s^2
//  pbar /= (rhoinfty * pow(cinfty, 2));
  double pbar = 1.0 / this->gamma_specificheat;

  // jet state
  double Tj = 600.0; // jet temperature, 600 K
  Tj /= ((this->gamma_specificheat - 1.0) * Tinfty); // non-dimensional jet temperature
  double cj = sqrt((this->gamma_specificheat - 1.0) * Tj); // non-dimensional jet speed of sound
  double Mj = 0.756;
  double uj = Mj * cj; // non-dimensional jet velocity
  double rhoj = pbar * this->gamma_specificheat / (this->gamma_specificheat - 1.0) / Tj;
  if (mpi::irank == 0) {

    std::cout << std::endl;
    std::cout << "Variables for the jet: " << std::endl;
    std::cout << "  Jet temperature (non-dimensional): " << Tj
              << "; ambient temperature (non-dimensional): " << 1.0 / (this->gamma_specificheat - 1.0) << std::endl;
    std::cout << "  Jet speed of sound (non-dimensional): " << cj
              << "; (gamma R T)^0.5 / cinfty: " << sqrt(1.4 * 287.0 * 600) / cinfty << std::endl;
    std::cout << "  Jet velocity (non-dimensional): " << uj << std::endl;
    std::cout << "  Jet density (non-dimensional): " << rhoj << std::endl;
    std::cout << std::endl;

  } // mpi::irank

  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_block = k - mygrid->iso[ZETA] + mygrid->iso_in_parent[ZETA];
    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_block = j - mygrid->iso[ETA] + mygrid->iso_in_parent[ETA];
      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_block = i - mygrid->iso[XI] + mygrid->iso_in_parent[XI];

        int l0 = mygrid->idx1D(i, j, k);

        // solution variables
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol[ivar])[l0] = 0.0;

        // base (or mean) state
        double ubar = uj * exp(-log(2.0) * pow(mygrid->cell[l0].xyz[YDIR], 2)); // b is the reference length
//        ubar = 0.0; // this gives an ambient base flow
//        ubar = uj; // this gives a uniform base flow
        (this->sol_mean[IVAR_UX])[l0] = ubar;
        (this->sol_mean[IVAR_UY])[l0] = 0.0; // a parallel jet
        (this->sol_mean[IVAR_UZ])[l0] = 0.0; // a two-dimensional jet
        (this->sol_mean[IVAR_P])[l0] = pbar; // uniform mean pressure

        // mean density
        double rhobar = -0.5 * (this->gamma_specificheat - 1.0) * (ubar - uj) * ubar
                      + ubar / (rhoj * uj) + (uj - ubar) / uj; // eq. (19) of Agarwal et al. (AIAA J. Vol. 42, No. 1, 2004)
        rhobar = 1.0 / rhobar;

        // mean entropy
        double sbar = (0.0 - this->gamma_specificheat * log(rhobar)) / this->gamma_specificheat; // pressure is uniform (pbar = pinfty)
        //
        (this->sol_mean[IVAR_S])[l0] = sbar;

      } // i
    } // j
  } // k

  return;

} // State::initialize_state_linearizedEuler_2Djet_with_a_harmonic_source



void State::initialize_state_linearizedEuler_cylinderScattering(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // harmonic, spatially distributed acoustic source emitting sound which is 
  // scattered by multiple rigid cylinders
  // see the fourth computational aeroacoustics workshop on benchmark problems (2004)
  // also see Sherer (J. Acoust. Soc. Am. 2004)

  double pbar = 1.0 / this->gamma_specificheat;
  double rhobar = 1.0;
  double Tbar = 1.0 / (this->gamma_specificheat - 1.0);

  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_block = k - mygrid->iso[ZETA] + mygrid->iso_in_parent[ZETA];

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_block = j - mygrid->iso[ETA] + mygrid->iso_in_parent[ETA];

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_block = i - mygrid->iso[XI] + mygrid->iso_in_parent[XI];

        int l0 = mygrid->idx1D(i, j, k);

        // solution variables
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol[ivar])[l0] = 0.0;

        // base (or mean) state
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol_mean[ivar])[l0] = 0.0;
        (this->sol_mean[IVAR_P])[l0] = pbar;

      } // i
    } // j
  } // k

  return;

} // State::initialize_state_linearizedEuler_cylinderScattering



void State::initialize_state_linearizedEuler_TannaTPN49(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  double pbar = 1.0 / this->gamma_specificheat;
  double rhobar = 1.0;
  double Tbar = 1.0 / (this->gamma_specificheat - 1.0);

  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_block = k - mygrid->iso[ZETA] + mygrid->iso_in_parent[ZETA];

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_block = j - mygrid->iso[ETA] + mygrid->iso_in_parent[ETA];

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_block = i - mygrid->iso[XI] + mygrid->iso_in_parent[XI];

        int l0 = mygrid->idx1D(i, j, k);

        // solution variables
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol[ivar])[l0] = 0.0;

        // base (or mean) state
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol_mean[ivar])[l0] = 0.0;
        (this->sol_mean[IVAR_P])[l0] = pbar;

      } // i
    } // j
  } // k

  return;

} // State::initialize_state_linearizedEuler_TannaTPN49



void State::initialize_state_linearizedEuler_KBKCombustor(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  double pbar = 1.0 / this->gamma_specificheat;
  double rhobar = 1.0;
  double Tbar = 1.0 / (this->gamma_specificheat - 1.0);

  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_block = k - mygrid->iso[ZETA] + mygrid->iso_in_parent[ZETA];

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_block = j - mygrid->iso[ETA] + mygrid->iso_in_parent[ETA];

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_block = i - mygrid->iso[XI] + mygrid->iso_in_parent[XI];

        int l0 = mygrid->idx1D(i, j, k);

        // solution variables
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol[ivar])[l0] = 0.0;

        // base (or mean) state
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol_mean[ivar])[l0] = 0.0;
        (this->sol_mean[IVAR_P])[l0] = pbar;

      } // i
    } // j
  } // k

  return;

} // State::initialize_state_linearizedEuler_KBKCombustor



void State::initialize_state_linearizedEuler_linearNozzle(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  double pbar = 1.0 / this->gamma_specificheat;
  double rhobar = 1.0;
  double Tbar = 1.0 / (this->gamma_specificheat - 1.0);

  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_block = k - mygrid->iso[ZETA] + mygrid->iso_in_parent[ZETA];

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_block = j - mygrid->iso[ETA] + mygrid->iso_in_parent[ETA];

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_block = i - mygrid->iso[XI] + mygrid->iso_in_parent[XI];

        int l0 = mygrid->idx1D(i, j, k);

        // solution variables
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol[ivar])[l0] = 0.0;

        // base (or mean) state
        for (int ivar = IVAR_S; ivar <= IVAR_P; ivar++)
          (this->sol_mean[ivar])[l0] = 0.0;
        (this->sol_mean[IVAR_P])[l0] = pbar;

      } // i
    } // j
  } // k

  if (myinput->num_scalar > 0) {
    int ivar_shift = IVAR_P + 1;

    for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
      int k_in_block = k - mygrid->iso[ZETA] + mygrid->iso_in_parent[ZETA];

      for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
        int j_in_block = j - mygrid->iso[ETA] + mygrid->iso_in_parent[ETA];

        for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
          int i_in_block = i - mygrid->iso[XI] + mygrid->iso_in_parent[XI];

          int l0 = mygrid->idx1D(i, j, k);

          for (int ivar = 0; ivar < myinput->num_scalar; ivar++) {
            (this->sol[ivar_shift+ivar])[l0] = 0.0;
            (this->sol_mean[ivar_shift+ivar])[l0] = 1.0;
          } // ivar

          double exp_factor = -log(4.0) / (0.1 * 0.1);
          double x0 = 0.1, y0 = 0.0;
          double xs = mygrid->cell[l0].xyz[XDIR] - x0;
          double ys = mygrid->cell[l0].xyz[YDIR] - y0;
          double amplitude = 0.01;
          (this->sol[IVAR_P+1])[l0] = amplitude * exp(exp_factor * (xs*xs + ys*ys));

        } // i
      } // j
    } // k
  } // myinput->num_scalar

  return;

} // State::initialize_state_linearizedEuler_linearNozzle



void State::backup_the_current_solution() {

  for (int ivar = 0; ivar < this->num_vars_sol; ivar++) {
    for (int l0 = 0; l0 < this->num_samples; l0++) {

      (this->sol_old[ivar])[l0] = (this->sol[ivar])[l0];

    } // l0
  } // ivar

  return;

} // State::backup_the_current_solution



void State::compute_dependent_variables(double **sol_cur) {

  // physical model
  if (this->model_pde == "LINEAR_ACOUSTICS")
    this->compute_dependent_variables_acoustics(sol_cur);

  else if (this->model_pde == "LEE" ||
           this->model_pde == "LEE_SCALAR") // auxiliary variables are the same as those for LEE
    this->compute_auxiliary_variables_linearizedEuler(sol_cur);

  else if (this->model_pde == "LEE_MIXFRAC_CONSTGAMMA")
    this->compute_auxiliary_variables_linearizedEuler_mixfrac_constgamma(sol_cur);

  else
    mpi::graceful_exit("This is a simulation for a unknown physical model.");

  return;

} // compute_dependent_variables



void State::compute_dependent_variables_acoustics(double **sol_cur) {

  double c_0 = 1.0; // ambient speed of sound

  // density fluctuation (from isentropic speed of sound relation)
  for (int l0 = 0; l0 < this->num_samples; l0++)
    (this->sol[IVAR_RHO])[l0] = (sol_cur[IVAR_P])[l0] / pow(c_0, 2.0);

  return;

} // compute_dependent_variables_acoustics



void State::compute_auxiliary_variables_linearizedEuler(double **sol_cur) {

  double gamma = this->gamma_specificheat;
  double gammaInv = 1.0 / gamma;
  double gammaOverGammaMinus1 = gamma / (gamma - 1.0);

  for (int l0 = 0; l0 < this->num_samples; l0++) {

    double sbar = (this->sol_mean[IVAR_S])[l0];
    double pbar = (this->sol_mean[IVAR_P])[l0];
    double pbarInv = 1.0 / pbar;
    double rhobar = pow( gamma * pbar, gammaInv) * exp(-sbar);
    double rhobarInv = 1.0 / rhobar;
    double Tbar = pbar * gammaOverGammaMinus1 * rhobarInv;

    double sPrime = (sol_cur[IVAR_S])[l0];
    double pPrime = (sol_cur[IVAR_P])[l0];

    // mean density
    (this->sol_aux[IAUX_RHO_MEAN])[l0] = rhobar;

    // density fluctuation (from linearized entropy expression)
    (this->sol_aux[IAUX_RHO])[l0] = gammaInv * (pPrime*pbarInv)
                                  - sPrime;
    (this->sol_aux[IAUX_RHO])[l0] *= rhobar;

    // mean temperature
    (this->sol_aux[IAUX_T_MEAN])[l0] = Tbar;

    // temperature fluctuation (from linearized equation of state)
    (this->sol_aux[IAUX_T])[l0] = (pPrime*pbarInv)
                                - (this->sol_aux[IAUX_RHO])[l0] * rhobarInv;
    (this->sol_aux[IAUX_T])[l0] *= Tbar;

  } // l0

  return;

} // compute_auxiliary_variables_linearizedEuler



void State::compute_auxiliary_variables_linearizedEuler_mixfrac_constgamma(double **sol_cur) {

  double gamma = this->gamma_specificheat;
  double gammaInv = 1.0 / gamma;
  double gammaMinus1OverGamma = (gamma - 1.0)/gamma;

  for (int l0 = 0; l0 < this->num_samples; l0++) {

    double sbar = (this->sol_mean[IVAR_S])[l0];
    double pbar = (this->sol_mean[IVAR_P])[l0];
    double pbarInv = 1.0 / pbar;
    double rhobar = (this->sol_aux[IAUX_RHO_MEAN])[l0];
    double rhobarInv = 1.0 / rhobar;
    double Tbar = (this->sol_aux[IAUX_T_MEAN])[l0];

    double sPrime = (sol_cur[IVAR_S])[l0];
    double pPrime = (sol_cur[IVAR_P])[l0];
    double ZPrime = (sol_cur[IVAR_Z])[l0];

    // density fluctuation (from linearized entropy expression)
    (this->sol_aux[IAUX_RHO])[l0] = gammaInv * (pPrime*pbarInv)
                                  - sPrime/(this->sol_aux[IAUX_CP])[l0]
                                  - (this->sol_aux[IAUX_PSI])[l0] * ZPrime;
    (this->sol_aux[IAUX_RHO])[l0] *= rhobar;

    // temperature fluctuation (from linearized equation of state)
    (this->sol_aux[IAUX_T])[l0] = sPrime/(this->sol_aux[IAUX_CP])[l0]
                                + gammaMinus1OverGamma * (pPrime*pbarInv)
                                - ((this->sol_aux[IAUX_DCPDZ])[l0]/(this->sol_aux[IAUX_CP])[l0] - (this->sol_aux[IAUX_PSI])[l0]) * ZPrime;
    (this->sol_aux[IAUX_T])[l0] *= Tbar;

  } // l0

  return;

} // compute_auxiliary_variables_linearizedEuler_mixfrac_constgamma



void State::prescribe_on_boundary_solution(Geometry::StructuredBoundaryCondition *myboundary, Geometry::StructuredGrid *mygrid, double time, double **myboundarydata) {

  // physical model
  if (this->model_pde == "LINEAR_ACOUSTICS")
    this->prescribe_on_boundary_solution_acoustics(myboundary, mygrid, time, myboundarydata);

  else if (this->model_pde == "LEE") {
    std::cout << "Dirichlet boundary for the linearized Euler model is not implemented yet." << std::endl;
    mpi::graceful_exit("Dirichlet boundary for the linearized Euler model is not implemented yet.");
  } // this->model_pde

  else if (this->model_pde == "LEE_SCALAR") {
    std::cout << "Dirichlet boundary for the linearized Euler model with passive scalars is not implemented yet." << std::endl;
    mpi::graceful_exit("Dirichlet boundary for the linearized Euler model with a scalar is not implemented yet.");
  } // this->model_pde

  else if (this->model_pde == "LEE_MIXFRAC_CONSTGAMMA") {
    std::cout << "Dirichlet boundary for the linearized Euler model with mixture fraction and a constant \\gamma is not implemented yet." << std::endl;
    mpi::graceful_exit("Dirichlet boundary for the linearized Euler model with a scalar is not implemented yet.");
  } // this->model_pde

  else
    mpi::graceful_exit("This is a simulation for a unknown physical model.");

  return;

} // State::prescribe_on_boundary_solution



void State::prescribe_on_boundary_solution_acoustics(Geometry::StructuredBoundaryCondition *myboundary, Geometry::StructuredGrid *mygrid, double time, double **myboundarydata) {

  if (this->simulation == "CASE_PLANE_WAVE")
    {
    int idir_wave_propagation = XDIR;

    for (int kb = myboundary->is[ZETA]; kb <= myboundary->ie[ZETA]; kb++) {
      int k_in_grid = kb - myboundary->is[ZETA] + myboundary->is_in_parent[ZETA];

      for (int jb = myboundary->is[ETA]; jb <= myboundary->ie[ETA]; jb++) {
        int j_in_grid = jb - myboundary->is[ETA] + myboundary->is_in_parent[ETA];

        for (int ib = myboundary->is[XI]; ib <= myboundary->ie[XI]; ib++) {
          int i_in_grid = ib - myboundary->is[XI] + myboundary->is_in_parent[XI];

          int lb = myboundary->idx1D(ib, jb, kb);
          int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

          if (mygrid->cell[l0].iblank != BLANKED) { // hole points are bypassed

            double *plane_wave = this->acoustics_plane_wave(time, mygrid->cell[l0].xyz[idir_wave_propagation], idir_wave_propagation);

            for (int ivar = IVAR_RHO; ivar <= IVAR_P; ivar++)
              (myboundarydata[ivar])[lb] = plane_wave[ivar];

            DEALLOCATE_1DPTR(plane_wave);

          } // mygrid->cell[l0].iblank
        } // ib
      } // jb
    } // kb
    }
  else
    mpi::graceful_exit("Unknown simulation for the current physical model.");

  return;

} // State::prescribe_on_boundary_solution_acoustics
