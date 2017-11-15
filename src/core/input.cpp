#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "input.h"

UserInput::UserInput() {


} // UserInput::UserInput



UserInput::~UserInput(void) {

} // UserInput::~UserInput



void UserInput::set(int argc, char * argv[]) {

  std::string str_stdout;
  std::stringstream str_counter;

  // read input-deck data and if needed, manually set something
  set_inputDeck(argc, argv);
  set_manual(argc, argv);



  // inspect the input data
  check_consistency_between_physical_model_and_simulation();

  // additional works using the input data
  get_number_of_variables();

  // spatial discretization
  if (OA_spatial > MAX_ORDER_ACCURACY)
    mpi::graceful_exit("MAX_ORDER_ACCURACY is exceeded for finite difference.");
  get_number_of_ghostCells_due_finiteDifference();

  // low-pass filter
  //if (OA_filter =< OA_spatial)
  //  mpi::graceful_exit("It is generally recommended to use a filter whose order of accuracy is higher than the finite difference.");
  if (filter_strength < 0.0 || filter_strength > 1.0)
    mpi::graceful_exit("The filter strength is out-of-bound.");
  if (filter_scheme == "STANDARD_CENTRAL") {
    MESSAGE_STDOUT("For standard central filter, filter strength is always the maximum; i.e. 1.0");
    filter_strength = 1.0;
  } // filter_scheme
  if (filter_blend < 0.0 || filter_blend > 1.0)
    mpi::graceful_exit("The filter-blending parameter is out-of-bound.");
  get_number_of_ghostCells_due_filter();

  // overset
  get_number_of_ghostCells_due_overset();

  // ghost cell
  num_cells_ghost = std::max( std::max(num_cells_ghost_finiteDifference, num_cells_ghost_filter), num_cells_ghost_overset);
  //
  str_counter.str("");
  str_counter << std::setw(6) << num_cells_ghost;
  str_stdout = "Number of ghost cells is " + str_counter.str();
  MESSAGE_STDOUT(str_stdout);

  // domain decomposition
  check_consistency_domainDecomposition();

  // files
  std::string file_base_local(file_base); // file_base is defined at core/param.h
  file_grid = file_base_local + ".xyz";
  file_metrics = file_base_local + ".met.q";
  file_solution = file_base_local; // contains only a base name; its extension is determined by the type of file
  file_varname_metrics = file_base_local + ".met.nam";
  file_varname_solution = file_base_local + ".nam";

  // determine whether we restart or not
  present_file_solution_in = FALSE;
  file_solution_in = "not-a-solution-file";
  if (argc >= 2) { // if there is at least one command-line argument, the first one is always a restart solution file

    present_file_solution_in = TRUE;
    file_solution_in = argv[1]; // argv[0] is reserved for the executable

  } // argc

  if (harmonicWave_waveType != "NONE")
    check_consistency_wave();

  return;

} // UserInput::set



void UserInput::set_inputDeck(int argc, char * argv[]) {

  std::string dummy;

  // geometry
  inputDeck::get_userInput("DIMENSION",num_dim);
  inputDeck::get_userInput("AXISYMMETRY",dummy); axisym = truefalse_2int(dummy);
  if (axisym == TRUE) {
    inputDeck::get_userInput("AXISYMMETRY","AXIS_OF_SYMMETRY",dummy); axis_of_sym = xyz_2int(dummy);
  } // axisym

  // physical model
  inputDeck::get_userInput("PHYSICAL_MODEL",model_pde);
  inputDeck::get_userInput("FLUID_MODEL",model_fluid);
  if (model_pde == "LEE_SCALAR") {
    inputDeck::get_userInput("PHYSICAL_MODEL","NUMBER_SCALAR",num_scalar);
    if (!(num_scalar >= 1))
      mpi::graceful_exit("For PHYSICAL_MODEL = " + model_pde + ", at least one scalar should be solved for.");
  } // model_pde
  else if (model_pde == "LEE_MIXFRAC_CONSTGAMMA")
    num_scalar = 1; // 1 for mixture fraction fluctuation Z'
  else
    num_scalar = 0;

  // simulation
  inputDeck::get_userInput("SIMULATION",simulation);

  // thermodynamics
  inputDeck::get_userInput("GAMMA_SPECIFICHEAT",gamma_specificheat);
  assert(gamma_specificheat > 0.0);
  inputDeck::get_userInput("SPECIFICHEAT_P",c_p);
  assert(c_p > 0.0);

  // temporal discretization
  inputDeck::get_userInput("FIX_DT_OR_CFL",fix_dt_or_cfl);
  dt = 1E-6; // a dummy time-step size
  if (fix_dt_or_cfl == "FIX_DT")
    inputDeck::get_userInput("DT",dt);
  else if (fix_dt_or_cfl == "FIX_CFL")
    inputDeck::get_userInput("CFL",cfl);
  else
    mpi::graceful_exit("You may fix either DT or CFL in time integration.");
  inputDeck::get_userInput("TEMPORAL_SCHEME",temporal_scheme);
  inputDeck::get_userInput("NUM_TOTAL_TIMESTEP",num_time_steps);
  inputDeck::get_userInput("REPORT_TIMESTEP",report_freq);

  // spatial discretization
  inputDeck::get_userInput("SPATIAL_SCHEME",spatial_scheme);
  if (spatial_scheme == "STANDARD_CENTRAL")
    inputDeck::get_userInput("SPATIAL_SCHEME","ORDER_ACCURACY",OA_spatial);

  // low-pass filter
  inputDeck::get_userInput("SPATIAL_FILTER",filter_scheme);
  if (filter_scheme == "NONE")
    do_filter = FALSE;
  else
    do_filter = TRUE;
  if (do_filter == TRUE) {
    if (filter_scheme == "STANDARD_CENTRAL")
      inputDeck::get_userInput("SPATIAL_FILTER","ORDER_ACCURACY",OA_filter);

    else if (filter_scheme == "SFO11P")
      inputDeck::get_userInput("SPATIAL_FILTER","STRENGTH",filter_strength);

    else
      mpi::graceful_exit("The filter " + filter_scheme + " is not supported.");

    inputDeck::get_userInput("SPATIAL_FILTER","BLEND",filter_blend);
  } // do_filter

  // buffer zone
  inputDeck::get_userInput("BUFFER_ZONE","POLYNOMIAL_ORDER",buffer_polynomial_order);
  inputDeck::get_userInput("BUFFER_ZONE","CONSTANT",buffer_constant);

  // metrics
  inputDeck::get_userInput("METRIC_SCHEME",scheme_metrics);

  // overset
  inputDeck::get_userInput("OVERSET",overset_format);
  if (overset_format == "NONE")
    do_overset = FALSE;
  else
    do_overset = TRUE;
  //
  if (do_overset == TRUE) {
    if (overset_format == "OVERTURE") {

      inputDeck::get_userInput("OVERSET","TYPE",overset_type_of_interpolation);
      inputDeck::get_userInput("OVERSET","ORDER_ACCURACY",overset_accuracy);

    } // overset_format
    else
      mpi::graceful_exit(overset_format+" is a unsupported format for overset-grid interpolation.");
  } // do_overset

  // domain decomposition
  inputDeck::get_userInput("DOMAIN_DECOMPOSITION",how2decompDomain);
  if (how2decompDomain == "DECOMP_FROMFILE")
    inputDeck::get_userInput("DOMAIN_DECOMPOSITION","FILE",file_decompDomain);

  // files
  inputDeck::get_userInput("FILE_FORMAT",type_file);
  //
  present_file_grid_in = TRUE;
  inputDeck::get_userInput("GRID_FILE",file_grid_in);

  if (do_overset == TRUE)
    inputDeck::get_userInput("OVERSET_FILE",file_overset);

  present_file_mean_in = FALSE;
  if (model_pde == "LEE" ||
      model_pde == "LEE_SCALAR" ||
      model_pde == "LEE_MIXFRAC_CONSTGAMMA" ||
      model_pde == "LNS") {

    present_file_mean_in = TRUE; // to solve linearized equations, mean (or base) states are needed
    inputDeck::get_userInput("BASESTATE_FILE",file_mean_in);

  } // model_pde
  present_file_aux_in = FALSE;
  if (model_pde == "LEE_MIXFRAC_CONSTGAMMA") {

    present_file_aux_in = TRUE; // some problems require to specify additional quantities used for simulation
                                // put them in auxiliary variables
    inputDeck::get_userInput("AUXVAR_FILE",file_aux_in);

  } // model_pde
  inputDeck::get_userInput("BC_FILE",file_boundary);

  // solution writing
  inputDeck::get_userInput("SOLUTION_WRITING_TIME",time_writing_solutions);

  // data probing
  do_probe == FALSE;
  num_probes = inputDeck::count_inputDeck_name("PROBE");
  if (num_probes > 0) {
    do_probe = TRUE;

    tmp_probe_name = new std::string[num_probes];
    ALLOCATE1D_INT_1ARG(tmp_probe_interval,num_probes); 
    ALLOCATE2D_DOUBLE(tmp_probe_xyz,num_probes, num_dim);
    double *xyz;
    ALLOCATE1D_DOUBLE_1ARG(xyz,DIM_MAX);

    // temporarily store probe-related data contained in the input deck so that 
    // they can be processed later
    for (int iprobe = 0; iprobe < num_probes; iprobe++) {
      inputDeck::get_userInput("PROBE","NAME",tmp_probe_name[iprobe],iprobe);
      inputDeck::get_userInput("PROBE","INTERVAL",tmp_probe_interval[iprobe],iprobe);
      inputDeck::get_userInput("PROBE","XYZ",DIM_MAX,xyz,iprobe);
      for (int idir = XDIR; idir < DIM_MAX; idir++)
        tmp_probe_xyz[iprobe][idir] = xyz[idir];
    } // iprobe
    DEALLOCATE_1DPTR(xyz);
  } // num_probes

  // time-harmonic wave parameters, if used
  inputDeck::get_userInput("HARMONIC_WAVE",harmonicWave_waveType);
  if (harmonicWave_waveType != "NONE") {

    if (harmonicWave_waveType == "WAVE_ACOUSTIC") {

      inputDeck::get_userInput("HARMONIC_WAVE","AMPLITUDE_RELATIVE",harmonicWave_amplitude);
      inputDeck::get_userInput("HARMONIC_WAVE","SHAPE",harmonicWave_waveForm);
      inputDeck::get_userInput("HARMONIC_WAVE","DIRECTION",dummy); harmonicWave_idir_propagation = xyz_2int(dummy);

    } // harmonicWave_waveType
    else if (harmonicWave_waveType == "WAVE_PRESSURE") {

      inputDeck::get_userInput("HARMONIC_WAVE","AMPLITUDE",harmonicWave_amplitude);
      inputDeck::get_userInput("HARMONIC_WAVE","SHAPE",harmonicWave_waveForm);
      inputDeck::get_userInput("HARMONIC_WAVE","DIRECTION",dummy); harmonicWave_idir_propagation = xyz_2int(dummy);

    } // harmonicWave_waveType
    else if (harmonicWave_waveType == "WAVE_ENTROPY") {

      inputDeck::get_userInput("HARMONIC_WAVE","AMPLITUDE",harmonicWave_amplitude);
      inputDeck::get_userInput("HARMONIC_WAVE","SHAPE",harmonicWave_waveForm);
      inputDeck::get_userInput("HARMONIC_WAVE","DIRECTION",dummy); harmonicWave_idir_propagation = xyz_2int(dummy);

    } // harmonicWave_waveType
    else if (harmonicWave_waveType == "WAVE_MIXFRAC") {

      inputDeck::get_userInput("HARMONIC_WAVE","AMPLITUDE",harmonicWave_amplitude);
      inputDeck::get_userInput("HARMONIC_WAVE","SHAPE",harmonicWave_waveForm);
      inputDeck::get_userInput("HARMONIC_WAVE","DIRECTION",dummy); harmonicWave_idir_propagation = xyz_2int(dummy);

    } // harmonicWave_waveType
    else
      mpi::graceful_exit("HARMONIC_WAVE = " + harmonicWave_waveType + " is not supported.");

    if (harmonicWave_waveForm == "WAVEFORM_PLANE") {

      if (harmonicWave_waveType != "WAVE_ACOUSTIC")
        mpi::graceful_exit("Currently, plane wave should be acoustic; use WAVEFORM_PLANE only with WAVE_ACOUSTIC.");
      inputDeck::get_userInput("HARMONIC_WAVE","WAVELENGTH",harmonicWave_wavelength);

    } // harmonicWave_waveForm
    else if (harmonicWave_waveForm == "WAVEFORM_HOMOGENEOUS") {

      inputDeck::get_userInput("HARMONIC_WAVE","PERIOD",harmonicWave_period);

    } // harmonicWave_waveForm
    else if (harmonicWave_waveForm == "WAVEFORM_GAUSSIAN_HALFWIDTH") {

      inputDeck::get_userInput("HARMONIC_WAVE","PERIOD",harmonicWave_period);
      inputDeck::get_userInput("HARMONIC_WAVE","HALFWIDTH",harmonicWave_halfWidth);

    } // harmonicWave_waveForm
    else if (harmonicWave_waveForm == "WAVEFORM_CUSTOM") {

      inputDeck::get_userInput("HARMONIC_WAVE","PERIOD",harmonicWave_period);

    } // harmonicWave_waveForm
  } // harmonicWave_waveType

  inputDeck::get_userInput("INFLOW_EXTERNAL",inflow_external);
  if (inflow_external != "NONE") {
    if (inflow_external == "TEMPORAL") { // only temporal variation comes from a file
                                         // thus, spatial information is prescribed inside the code
      inputDeck::get_userInput("INFLOW_EXTERNAL","INFLOW_FILE",file_inflow);
      inputDeck::get_userInput("INFLOW_EXTERNAL","ORDER_ACCURACY_INTERP_TIME",OA_time_inflow);
      inputDeck::get_userInput("INFLOW_EXTERNAL","SHAPE_SPACE",shape_space_inflow);

    } // inflow_external
    else
      mpi::graceful_exit(inflow_external + " is not supported for external-inflow option.");

  } // inflow_external
std::cout << inflow_external << std::endl;
std::cout << file_inflow << std::endl;
std::cout << OA_time_inflow << std::endl;
std::cout << shape_space_inflow << std::endl;
 mpi::graceful_exit("Bye bye!");

  // solution interpolation
  inputDeck::get_userInput("INTERPOLATE_SOLUTION",interp_fromWhichFormat);
  todos::add("Checking INTERPOLATE_SOLUTION should be bypassed for the other runs.");
  todos::add("Fix why overset and filtering crash while interpolating solutions.");
  if (interp_fromWhichFormat != "NONE") {

    // somehow, solution interpolation does not work when both overset and filtering are enabled at the same time
    // since filtering is necessary for solution interpolation, turn off overset capability for now (this has to be fixed)
    do_overset = FALSE;
    MESSAGE_STDOUT("Solution interpolation requires to disable overset.");

    present_file_mean_in = FALSE;
    present_file_aux_in = FALSE;
    MESSAGE_STDOUT("While interpolating solutions, stop reading a base-state or auxiliary file for linearized analysis.");

    // data type by which source data are written
    if (interp_fromWhichFormat == "TECPLOT_FE") {

      inputDeck::get_userInput("INTERPOLATE_SOLUTION","NUM_DIM_SOURCE",num_dim_interpSource); // e.g. 2 if the source data are 2-D
      inputDeck::get_userInput("INTERPOLATE_SOLUTION","NUM_VARS_SOURCE",num_vars_interpSource); // number of solution variables to be interpolated
      inputDeck::get_userInput("INTERPOLATE_SOLUTION","VARS_SOURCE",vars_interpSource);
      inputDeck::get_userInput("INTERPOLATE_SOLUTION","NUM_ZONES_SOURCE",num_zones_interpSource); // in case the source data consist of multiple zones with different element types

      if (num_dim_interpSource != 2)
        mpi::graceful_exit("NUM_DIM_SOURCE in INTERPOLATE_SOLUTION has been tested only for 2-D source data.");

    } // interp_fromWhichFormat
    else
      mpi::graceful_exit("Data format " + interp_fromWhichFormat + " is not supported for solution interpolation.");

    // reference scales for non-dimensionalization
    inputDeck::get_userInput("INTERPOLATE_SOLUTION","REF_LENGTH",scale_xyz);
    if (vars_interpSource == "PRHOU") { // p, \rho, u_i

      if (model_pde != "LEE")
        mpi::graceful_exit("PHYSICAL_MODEL " + model_pde + " is incompatible with variables " + vars_interpSource + " for interpolation.");

      inputDeck::get_userInput("INTERPOLATE_SOLUTION","REF_PRESSURE",scale_p);
      inputDeck::get_userInput("INTERPOLATE_SOLUTION","REF_DENSITY",scale_rho);
      inputDeck::get_userInput("INTERPOLATE_SOLUTION","REF_VELOCITY",scale_u);

    } // vars_interpSource
    else if (vars_interpSource == "PUS") { // p, u_i, s

      if (model_pde != "LEE")
        mpi::graceful_exit("PHYSICAL_MODEL " + model_pde + " is incompatible with variables " + vars_interpSource + " for interpolation.");

      inputDeck::get_userInput("INTERPOLATE_SOLUTION","REF_PRESSURE",scale_p);
      inputDeck::get_userInput("INTERPOLATE_SOLUTION","REF_VELOCITY",scale_u);
      inputDeck::get_userInput("INTERPOLATE_SOLUTION","REF_ENTROPY",scale_s);

    } // vars_interpSource
    else
      mpi::graceful_exit("VARS_SOURCE " + vars_interpSource + " is not supported for interpolation.");

    inputDeck::get_userInput("INTERPOLATE_SOLUTION","NUM_FILTER",num_filters_interpolation);

    inputDeck::get_userInput("INTERPOLATE_SOLUTION","WRITE_FORMAT",interpolate_into_which_PLOT3D);

    inputDeck::get_userInput("INTERPOLATE_SOLUTION","SOURCE_FILE",file_tecplot_ASCII);

  } // interp_fromWhichFormat

  // miscellaneous
  todos::add("Change the code so that dimensional reference quantities can be input.");
  todos::add("Put an option writing dt, cfl (total, convective, & viscous) in a plot3d file.");
  todos::add("Write a post-processor reading and averaging solution and function files.");
  todos::add("AUX in the zone header causes a trouble in reading Tecplot files.");
  todos::add("While interpolating solution, overset and filtering cannot be used simultaneously.");

  return;

} // UserInput::set_inputDeck



void UserInput::set_manual(int argc, char * argv[]) {

//  std::string str_stdout;
//  std::stringstream str_counter;

//  // geometry
//  num_dim = 2;
//  axisym = TRUE;
//  axis_of_sym = XDIR;

//  // physical model
//  model_pde = "LEE";
////  model_pde = "LINEAR_ACOUSTICS";
//  model_fluid = "IDEAL_GAS_CALORIC";

//  // simulation
//  simulation = "CASE_KBK_COMBUSTOR";
//  simulation = "CASE_TANNA_TPN49";
//  simulation = "CASE_SCATTERING_TWOCYLINDER";
//  simulation = "CASE_2DJET_WITH_A_HARMONIC_SOURCE";
//  simulation = "CASE_GAUSSIAN_PULSE";
//  simulation = "CASE_PLANE_WAVE";
//  check_consistency_between_physical_model_and_simulation();

//  // thermodynamics
//  gamma_specificheat = 1.4;

//  // variables
//  get_number_of_variables();

//  // temporal discretization
//  fix_dt_or_cfl = "FIX_DT"; // "FIX_CFL";
//  dt = 5.0E-3; //CASE_KBK_COMBUSTOR
//  dt = 8.0E-4; //CASE_TANNA_TPN49 (RANS jet base state)
//  dt = 2.0E-4; //CASE_TANNA_TPN49 (quiescent base state)
//  dt = 2.0E-3; //CASE_SCATTERING_TWOCYLINDER
//  dt = 1.0E-2; //CASE_2DJET_WITH_A_HARMONIC_SOURCE
//  dt = 2.0E-4; //CASE_PLANE_WAVE
//  cfl = 0.5;
//  temporal_scheme = "RUNGE_KUTTA4";
//  num_time_steps = 100000000; //CASE_KBK_COMBUSTOR
//  num_time_steps = 1250000; //CASE_TANNA_TPN49
//  num_time_steps = 20000; //CASE_SCATTERING_TWOCYLINDER
//  num_time_steps = 44160; //CASE_2DJET_WITH_A_HARMONIC_SOURCE
//  num_time_steps = 40000; //CASE_PLANE_WAVE
//  report_freq = 100; // report every report_freq time step

//  // finite difference
//  spatial_scheme = "STANDARD_CENTRAL";
//  spatial_scheme = "FDO11P";
//  OA_spatial = 4;
//  if (OA_spatial > MAX_ORDER_ACCURACY)
//    mpi::graceful_exit("MAX_ORDER_ACCURACY is exceeded for finite difference.");
//  get_number_of_ghostCells_due_finiteDifference();

//  // low-pass filter
//  do_filter = TRUE;
//  if (do_filter == TRUE) {

//    filter_scheme = "STANDARD_CENTRAL";
//    OA_filter = 8; //4;
//    filter_scheme = "SFO11P";
//    OA_filter = 2;
//    if (OA_filter =< OA_spatial)
//      mpi::graceful_exit("It is generally recommended to use a filter whose order of accuracy is higher than the finite difference.");

//    filter_strength = 0.2; //CASE_KBK_COMBUSTOR
//    //filter_strength = 0.6; //1.0; //CASE_TANNA_TPN49
//    if (filter_strength < 0.0 || filter_strength > 1.0)
//      mpi::graceful_exit("The filter strength is out-of-bound.");
//    if (filter_scheme == "STANDARD_CENTRAL") {
//      MESSAGE_STDOUT("For standard central filter, filter strength is always the maximum; i.e. 1.0");
//      filter_strength = 1.0;
//    } // filter_scheme
//
//    filter_blend = 1.0; //CASE_KBK_COMBUSTOR
//    //filter_blend = 1.0; //0.5; //CASE_TANNA_TPN49
//    if (filter_blend < 0.0 || filter_blend > 1.0)
//      mpi::graceful_exit("The filter-blending parameter is out-of-bound.");

//  } // do_filter
//  get_number_of_ghostCells_due_filter();

//  // buffer zone
//  buffer_polynomial_order = 2;
////  buffer_constant = 0.5;
//  buffer_constant = 5.0; //CASE_KBK_COMBUSTOR
////  buffer_constant = 10.0; //CASE_TANNA_TPN49 (RANS jet base state)
////  buffer_constant = 5.0; //CASE_TANNA_TPN49 (quiescent base state)
////  buffer_constant = 50.0; //CASE_PLANE_WAVE

//  // metrics
//  scheme_metrics = "THOMAS_LOMBARD";

//  // overset
//  do_overset = FALSE; //TRUE;
////  overset_format = "OVERTURE_NOTSCALING";
//  overset_format = "OVERTURE";
//  overset_type_of_interpolation = "EXPLICIT";
//  overset_accuracy = 4; //6;
//  get_number_of_ghostCells_due_overset();

//  // ghost cell
//  num_cells_ghost = std::max( std::max(num_cells_ghost_finiteDifference, num_cells_ghost_filter), num_cells_ghost_overset);
//  //
//  str_counter.str("");
//  str_counter << std::setw(6) << num_cells_ghost;
//  str_stdout = "Number of ghost cells is " + str_counter.str();
//  MESSAGE_STDOUT(str_stdout);

//  // domain decomposition
//  how2decompDomain = "DECOMP_1D";
////  how2decompDomain = "DECOMP_FROMFILE";
////  file_decompDomain = "/home/jeokim/Mesh/ASMENozzle/decomp.map"; //CASE_TANNA_TPN49
//  file_decompDomain = "/home/jeokim/Mesh/CTRSP16_Combustor/decomp.map"; //CASE_KBK_COMBUSTOR
//  check_consistency_domainDecomposition();

//  // files
//  type_file = "PLOT3D";
//  std::string file_base_local(file_base); // file_base is defined at core/param.h
//  file_grid = file_base_local + ".xyz";
//  file_metrics = file_base_local + ".met.q";
//  file_solution = file_base_local; // contains only a base name; its extension is determined by the type of file
//  file_varname_metrics = file_base_local + ".met.nam";
//  file_varname_solution = file_base_local + ".nam";

//  present_file_grid_in = TRUE;
//  file_grid_in = "/home/jeokim/Mesh/CTRSP16_Combustor/nozzle_binary.xyz";
////  file_grid_in = "/home/jeokim/Mesh/ASMENozzle/TannaNozzle_binary.xyz";
////  file_grid_in = "/home/jeokim/Mesh/cylinderScattering/twoCylinder_binary.xyz";
////  file_grid_in = "/home/jeokim/Mesh/jet_with_a_harmonic_source/jet2D_binary.xyz";
////  file_grid_in = "/home/jeokim/Mesh/oversetTests/twoBlock_binary.xyz";
//  if (present_file_grid_in == TRUE) {
//    MESSAGE_STDOUT("Grid file is " + file_grid_in + ".");
//  } // present_file_grid_in
//  else
//    mpi::graceful_exit("Specify a relevant grid file.");

// file_overset = "/home/jeokim/Mesh/CTRSP16_Combustor/overset.bin";
////  file_overset = "/home/jeokim/Mesh/ASMENozzle/overset.bin";
////  file_overset = "/home/jeokim/Mesh/cylinderScattering/overset_twoCylinder.bin";
////  file_overset = "/home/jeokim/Mesh/oversetTests/overset_twoBlock.bin";
//  if (do_overset == TRUE)
//    MESSAGE_STDOUT("Overset file is " + file_overset + ".");

//  // determine whether we restart or not
//  present_file_solution_in = FALSE;
//  file_solution_in = "not-a-solution-file";
//  if (argc >= 2) { // if there is at least one command-line argument, the first one is always a restart solution file
//
//    present_file_solution_in = TRUE;
//    file_solution_in = argv[1]; // argv[0] is reserved for the executable
//
//  } // argc
//  if ( present_file_solution_in == TRUE)
//    str_stdout = "Simulation restarts from " + file_solution_in + ".";
//  else
//    str_stdout = "Simulation starts from scratch.";
//  MESSAGE_STDOUT(str_stdout);

//  present_file_mean_in = FALSE; //TRUE; //CASE_KBK_COMBUSTOR
////  present_file_mean_in = TRUE; //CASE_TANNA_TPN49 (RANS jet base state)
////  present_file_mean_in = FALSE; //CASE_TANNA_TPN49 (quiescent base state)
////  file_mean_in = "/home/jeokim/work/TPN49/LEE/tecplotFE2plot3d/03.8thStdFilter_s1.0alpha1.0_100times/mean_camille.q"; //CASE_TANNA_TPN49 (RANS jet base state)
//  file_mean_in = "/home/jeokim/work/CTRSP16_Combustor/mean_from_LES/04.8thStdFilter_s1.0alpha1.0_100times/mean_camille.q"; //CASE_KBK_COMBUSTOR
//  if (present_file_mean_in == TRUE)
//    MESSAGE_STDOUT("Base-state file is " + file_mean_in + ".");

//  file_boundary = "/home/jeokim/Mesh/CTRSP16_Combustor/bc.in";
////  file_boundary = "/home/jeokim/Mesh/ASMENozzle/bc.in";
////  file_boundary = "/home/jeokim/Mesh/cylinderScattering/bc_twoCylinder.in";
////  file_boundary = "/home/jeokim/Mesh/jet_with_a_harmonic_source/bc.in";
////  file_boundary = "/home/jeokim/Mesh/oversetTests/bc.in";
//  MESSAGE_STDOUT("Boundary condition file is " + file_boundary + ".");

//  // solution writing
//  time_writing_solutions = 1.0E+01; //CASE_KBK_COMBUSTOR
////  time_writing_solutions = 1.0E+00; //CASE_TANNA_TPN49 (StD = 0.04)
////  time_writing_solutions = 1.0E-01; //CASE_TANNA_TPN49 (StD = 0.4)
////  time_writing_solutions = 2.0E-02; //CASE_SCATTERING_TWOCYLINDER
////  time_writing_solutions = 2E+00; //CASE_2DJET_WITH_A_HARMONIC_SOURCE
////  time_writing_solutions = 2E-01; //CASE_PLANE_WAVE

//  // time-harmonic wave parameters, if used
//  harmonicWave_waveType = "WAVE_ACOUSTIC"; //CASE_KBK_COMBUSTOR
//  harmonicWave_amplitude = 0.0001; //CASE_KBK_COMBUSTOR
//  harmonicWave_waveForm = "WAVEFORM_HOMOGENEOUS";
////  harmonicWave_waveType = "WAVE_ACOUSTIC"; //CASE_TANNA_TPN49
////  harmonicWave_amplitude = 0.00064517; //CASE_TANNA_TPN49 (0.06% of p_\infty; 130dB)
////  harmonicWave_waveForm = "WAVEFORM_PLANE";
////  harmonicWave_waveType = "WAVE_ENTROPY"; //CASE_TANNA_TPN49
////  harmonicWave_amplitude = 0.05; //CASE_TANNA_TPN49 (T^\prime/\bar{T} \approx s^\prime; 5% temperature fluctuation)
////  harmonicWave_waveForm = "WAVEFORM_HOMOGENEOUS";
//////  harmonicWave_waveForm = "WAVEFORM_GAUSSIAN_HALFWIDTH";
//  harmonicWave_halfWidth = 0.5;
//  harmonicWave_idir_propagation = XDIR;
//  harmonicWave_wavelength = 68.0; //CASE_KBK_COMBUSTOR
////  harmonicWave_wavelength = 70.0; //CASE_TANNA_TPN49 (Uj = 1.43): StD = 0.02
////  harmonicWave_wavelength = 35.0; //CASE_TANNA_TPN49 (Uj = 1.43): StD = 0.04
////  harmonicWave_wavelength = 3.5; //CASE_TANNA_TPN49 (Uj = 1.43): StD = 0.4 (jet column mode)
////  harmonicWave_wavelength = 1.9; //CASE_TANNA_TPN49 (Uj = 1.43): StD = 0.74 (wavelength slightly shorter the nozzle diameter)
////  harmonicWave_wavelength = 1.4; //CASE_TANNA_TPN49 (Uj = 1.43): StD = 1.0
////  harmonicWave_wavelength = 0.7; //CASE_TANNA_TPN49 (Uj = 1.43): StD = 2.0
//  check_consistency_wave();

//  // solution interpolation
//  //num_zones_interpSource = 2; // number of zones in the Tecplot file from which solutions are interpolated; count how many "ZONE T" in *.dat file //CASE_TANNA_TPN49
//  num_zones_interpSource = 4; // number of zones in the Tecplot file from which solutions are interpolated; count how many "ZONE T" in *.dat file //CASE_KBK_COMBUSTOR
//  num_dim_interpSource = 2; // 2 if 2-D
//  num_vars_interpSource = num_dim_interpSource + 2;
//  //vars_interpSource = "PRHOU"; //CASE_TANNA_TPN49
//  vars_interpSource = "PUS"; //CASE_KBK_COMBUSTOR
//  scale_xyz = 1.0;
//  scale_rho = 1.0 / 1.0792683;
//  scale_p = 1.0 / (1.0792683 * pow(356.59659, 2));
//  scale_u = 1.0 / 356.59659;
//  scale_s = 1.0 / 1005.0;
//  //
//  num_filters_interpolation = 0; //100;
//  interpolate_into_which_PLOT3D = "PLOT3D_FUNCTION";
//  //
//  //file_profile_FLUENT = "/home/jeokim/work/TPN49/RANS/03.UpTo72Dj/results_wRealizable_ke/pv.prof"; //CASE_TANNA_TPN49
//  //file_tecplot_ASCII = "/home/jeokim/work/TPN49/RANS/03.UpTo72Dj/results_wRealizable_ke/pv.dat"; //CASE_TANNA_TPN49
//  file_tecplot_ASCII = "/home/jeokim/work/CTRSP16_Combustor/computeMeanCombustorLES/04.average/avg_boundaryFixed.dat"; //CASE_KBK_COMBUSTOR

  return;

} // UserInput::set_manual



void UserInput::check_consistency_dimension() {

  // this routine checks a minimal consistency in the spatial dimension of the simulation

  if (axisym == TRUE) {

    if (num_dim != 2)
      mpi::graceful_exit("Axisymmetry requires 2-D.");

    if (axis_of_sym != XDIR)
      mpi::graceful_exit("The axis of symmetry should be x by definition of the coordinate system.");

  } // axisym

  return;

} // UserInput::check_consistency_dimension



void UserInput::check_consistency_between_physical_model_and_simulation() {

  // this routine checks a minimal consistency between the physical model and the simulation
  // make sure to add a relevant message below if a new physical model or a simulation is implemented

  int ok = NOT_OK;

  if (model_pde == "LINEAR_ACOUSTICS") {

    if ( simulation == "CASE_PLANE_WAVE" ) {

      MESSAGE_STDOUT("Linear acoustic equations are solved for plane-wave propagation.");
      ok = OK;

    } // simulation
    if ( simulation == "CASE_GAUSSIAN_PULSE" ) {

      MESSAGE_STDOUT("Linear acoustic equations are solved for an initially Gaussian pressure pulse.");
      ok = OK;

    } // simulation

  } // model_pde
  else if (model_pde == "LEE") {

    if ( simulation == "CASE_2DJET_WITH_A_HARMONIC_SOURCE" ) {

      MESSAGE_STDOUT("Linearized Euler equations are solved for acoustic wave refraction by a 2-D parallel jet.");
      ok = OK;

    } // simulation
    if ( simulation == "CASE_SCATTERING_TWOCYLINDER" ) {

      MESSAGE_STDOUT("Linearized Euler equations are solved for acoustic scattering by two rigid cylinders.");
      ok = OK;

    } // simulation
    if ( simulation == "CASE_TANNA_TPN49" ) {

      MESSAGE_STDOUT("Linearized Euler equations are solved for the subsonic hot jet of the Tanna's TPN 49.");
      ok = OK;

    } // simulation
    if ( simulation == "CASE_KBK_COMBUSTOR" ) {

      MESSAGE_STDOUT("Linearized Euler equations are solved for the KBK combustor set-up.");
      ok = OK;

    } // simulation
    if ( simulation == "CASE_LINEAR_NOZZLE" ) {

      MESSAGE_STDOUT("Linearized Euler equations are solved for the linear nozzle set-up.");
      ok = OK;

    } // simulation

  } // model_pde
  else if (model_pde == "LEE_SCALAR") {

    if ( simulation == "CASE_LINEAR_NOZZLE" ) {

      std::stringstream str_counter;
      str_counter << num_scalar;

      MESSAGE_STDOUT("Linearized Euler equations with " + str_counter.str() + " passive scalar(s) are solved for the linear nozzle set-up.");
      ok = OK;

    } // simulation

  } // model_pde
  else if (model_pde == "LEE_MIXFRAC_CONSTGAMMA") {

    if ( simulation == "CASE_LINEAR_NOZZLE" ) {

      MESSAGE_STDOUT("Linearized Euler equations with mixture fraction and a constant \\gamma are solved for the linear nozzle set-up.");
      ok = OK;

    } // simulation
    if ( simulation == "CASE_EXPONENTIAL_HORN" ) {

      MESSAGE_STDOUT("Linearized Euler equations with mixture fraction and a constant \\gamma are solved for the exponential horn set-up.");
      ok = OK;

    } // simulation

  } // model_pde
  else if (model_pde == "LNS") {

    if ( simulation == "CASE_KBK_COMBUSTOR" ) {

      MESSAGE_STDOUT("Linearized Navier--Stokes equations are solved for the KBK combustor set-up.");
      ok = OK;

    } // simulation

  } // model_pde

  if (!ok)
    mpi::graceful_exit("Unknown simulation for the current physical model.");

  return;

} // UserInput::check_consistency_between_physical_model_and_simulation



void UserInput::check_consistency_wave() {

  // this routine checks a minimal consistency in the time-harmonic wave 
  // specified as boundary condition or buffer zone

  if (harmonicWave_waveType == "NONE")
    return;

  if (harmonicWave_waveType == "WAVE_ACOUSTIC")
    if (harmonicWave_waveForm != "WAVEFORM_PLANE" &&
        harmonicWave_waveForm != "WAVEFORM_HOMOGENEOUS" &&
        harmonicWave_waveForm != "WAVEFORM_GAUSSIAN_HALFWIDTH" &&
        harmonicWave_waveForm != "WAVEFORM_CUSTOM")
      mpi::graceful_exit("SHAPE = " + harmonicWave_waveForm + " is a unknown wave form for HARMONIC_WAVE = " + harmonicWave_waveType);

  if (harmonicWave_waveType == "WAVE_PRESSURE")
    if (harmonicWave_waveForm != "WAVEFORM_HOMOGENEOUS" &&
        harmonicWave_waveForm != "WAVEFORM_GAUSSIAN_HALFWIDTH" &&
        harmonicWave_waveForm != "WAVEFORM_CUSTOM")
      mpi::graceful_exit("SHAPE = " + harmonicWave_waveForm + " is a unknown wave form for HARMONIC_WAVE = " + harmonicWave_waveType);

  else if (harmonicWave_waveType == "WAVE_ENTROPY")
    if (harmonicWave_waveForm != "WAVEFORM_HOMOGENEOUS" &&
        harmonicWave_waveForm != "WAVEFORM_GAUSSIAN_HALFWIDTH" &&
        harmonicWave_waveForm != "WAVEFORM_CUSTOM")
      mpi::graceful_exit("SHAPE = " + harmonicWave_waveForm + " is a unknown wave form for HARMONIC_WAVE = " + harmonicWave_waveType);

  else if (harmonicWave_waveType == "WAVE_MIXFRAC")
    if (harmonicWave_waveForm != "WAVEFORM_HOMOGENEOUS" &&
        harmonicWave_waveForm != "WAVEFORM_GAUSSIAN_HALFWIDTH" &&
        harmonicWave_waveForm != "WAVEFORM_CUSTOM")
      mpi::graceful_exit("SHAPE = " + harmonicWave_waveForm + " is a unknown wave form for HARMONIC_WAVE = " + harmonicWave_waveType);

  else
    mpi::graceful_exit("HARMONIC_WAVE = " + harmonicWave_waveType + " is not implemented.");

  if (harmonicWave_amplitude < 0.0)
    mpi::graceful_exit("The wave amplitude cannot be negative.");

  return;

} // UserInput::check_consistency_wave



void UserInput::check_consistency_domainDecomposition() {

  // this routine checks a minimal consistency in the domain-decomposition strategy

  std::ifstream ifs;

  if (how2decompDomain == "DECOMP_1D")
    mpi::wait_allothers("The domain (and each block) is decomposed along a single direction with the largest number of cells.");

  else if (how2decompDomain == "DECOMP_FROMFILE") {

    mpi::wait_allothers("Information on how to decompose the domain is read from a file.");

    // to avoid file-opening attempts by too many cores, only the global head core 
    // does the existence check
    if (mpi::irank == 0) {

      ifs.open(cstr_to_constchar(file_decompDomain), std::ifstream::in);
      if (!ifs.is_open())
        mpi::graceful_exit("The file for your domain decomposition does not exist.");
      ifs.close();

    } // mpi::irank
    mpi::wait_allothers();
  } // how2decompDomain
  else
    mpi::graceful_exit("The current domain-decomposition method is unknown.");

  return;

} // UserInput::check_consistency_wave



void UserInput::get_number_of_ghostCells_due_finiteDifference() {

  if (spatial_scheme == "STANDARD_CENTRAL")
    num_cells_ghost_finiteDifference = OA_spatial / 2;

  else if (spatial_scheme == "FDO11P")
    num_cells_ghost_finiteDifference = (FDO11P - 1) / 2;

  else
    mpi::graceful_exit("Unknown spatial scheme.");

  return;

} // UserInput::get_number_of_ghostCells_due_finiteDifference



void UserInput::get_number_of_ghostCells_due_filter() {

  num_cells_ghost_filter = 0;
  if (do_filter == FALSE)
    return;

  if (filter_scheme == "STANDARD_CENTRAL")
    num_cells_ghost_filter = OA_filter / 2;

  else if (filter_scheme == "SFO11P")
    num_cells_ghost_filter = (SFO11P - 1) / 2;

  else
    mpi::graceful_exit("Unknown filter scheme.");

  return;

} // UserInput::get_number_of_ghostCells_due_filter



void UserInput::get_number_of_ghostCells_due_overset() {

  num_cells_ghost_overset = 0;
  if (do_overset == FALSE)
    return;

  if (overset_type_of_interpolation == "EXPLICIT")
    num_cells_ghost_overset = overset_accuracy;

  else
    mpi::graceful_exit("Only explicit overset interpolation is supported for now.");

  return;

} // UserInput::get_number_of_ghostCells_due_overset



void UserInput::get_number_of_variables() {

  if (model_pde == "LINEAR_ACOUSTICS") {

    num_vars_sol = DIM_MAX + 2; // rho', u_i', p'
    num_vars_mean = 0;
    num_vars_meanGradient = 0;
    num_vars_aux = 0;
 
  } // model_pde
  else if (model_pde == "LEE") {

    num_vars_sol = DIM_MAX + 2; // s', u_i', p'
    num_vars_mean = num_vars_sol;
    num_vars_meanGradient = num_vars_mean;
    num_vars_aux = 2 * 2 + 1; // rho', T', and their means, c_p

  } // model_pde
  else if (model_pde == "LEE_SCALAR") {

    num_vars_sol = DIM_MAX + 2; // s', u_i', p'
    num_vars_sol += num_scalar; // passive scalar(s)
    num_vars_mean = num_vars_sol;
    num_vars_meanGradient = num_vars_mean;
    num_vars_aux = 2 * 2 + 1; // rho', T', and their means, c_p

  } // model_pde
  else if (model_pde == "LEE_MIXFRAC_CONSTGAMMA") {

    num_vars_sol = DIM_MAX + 2; // s', u_i', p'
    num_vars_sol += num_scalar; // mixture-fraction fluctuation Z'
    num_vars_mean = num_vars_sol;
    num_vars_meanGradient = num_vars_mean;
    num_vars_aux = 2 * 2 + 1; // rho', T', and their means, c_p
    num_vars_aux += 2; // dc_p/dZ, and \Psi

  } // model_pde
  else
    mpi::graceful_exit("PHYSICAL_MODEL = " + model_pde + " is unknown, and number of variables cannot be determined.");

  return;

} // UserInput::get_number_of_variables



int UserInput::truefalse_2int(std::string true_or_false) {

  str_toupper(true_or_false);

  int int2return = NONE;
  if (true_or_false == "TRUE")
    int2return = TRUE;
  else if (true_or_false == "FALSE")
    int2return = FALSE;
  else
    assert(0);

  return int2return;

} // UserInput::truefalse_2int



std::string UserInput::int_2truefalse(int true_or_false) {

  std::string str2return = "NONE";
  if (true_or_false == TRUE)
    str2return = "TRUE";
  else if (true_or_false = FALSE)
    str2return = "FALSE";
  else
    assert(0);

  return str2return;

} // UserInput::int_2truefalse



int UserInput::xyz_2int(std::string xyz) {

  str_toupper(xyz);

  int int2return = NONE;
  if (xyz == "X")
    int2return = XDIR;
  else if (xyz == "Y")
    int2return = YDIR;
  else if (xyz == "R")
    int2return = RDIR;
  else if (xyz == "Z")
    int2return = ZDIR;
  else
    assert(0);

  return int2return;

} // UserInput::xyz_2int



std::string UserInput::int_2xyz(int xyz) {

  std::string str2return = "NONE";
  if (xyz == XDIR)
    str2return = "x";
  else if (xyz == YDIR)
    str2return = "y";
  else if (xyz == ZDIR)
    str2return = "z";
  else
    assert(0);

  return str2return;

} // UserInput::int_2xyz



namespace inputDeck {

const char delimiter[] = " =,;:\t"; // delimiters used for input deck

int numLinesInputDeck = 0;
std::vector<std::string> linesInputDeck;

int numEntriesInputDeck = 0;
t_EntryInputDeck *entriesInputDeck;

void parse_linesInputDeck() {

  std::stringstream str_counter;

  // Step 1: if there is a comment symbol #, erase the remaining characters including # itself in the current line
  for (int iline = 0; iline < numLinesInputDeck; iline++) {
    for (int ichar = 0; ichar < linesInputDeck[iline].size(); ichar++) {
      if (linesInputDeck[iline].at(ichar) == '#') { // a comment symbol is detected
        linesInputDeck[iline].erase(linesInputDeck[iline].begin()+ichar,linesInputDeck[iline].end());
        break;

      } // linesInputDeck[iline].at(ichar)

    } // ichar
    //if (mpi::irank == 0)
    //  std::cout << "After comment symbols are removed: " << linesInputDeck[iline] << "<EOL>" << std::endl;
  } // iline



  // Step 2: if there is a backslash \, remove it, merge the current line into the next one, and empty the current line
  for (int iline = 0; iline < numLinesInputDeck; iline++) {
    int merge = FALSE; // by default, do not merge this line with the next one
    int ibackslash = NONE; // by default, there is no backslash anywhere

    for (int ichar = 0; ichar < linesInputDeck[iline].size(); ichar++) {
      if (linesInputDeck[iline].at(ichar) == '\\') { // a backslash is detected

        merge = TRUE;
        ibackslash = ichar;
        break;

      } // linesInputDeck[iline].at(ichar)
    } // ichar

    if (merge == TRUE) {
      assert(ibackslash != NONE);
      if (iline == numLinesInputDeck-1) { // if there is no more next line, just remove a backslash

        linesInputDeck[iline].erase(linesInputDeck[iline].begin()+ibackslash,linesInputDeck[iline].end());

      } // iline
      else {

        std::string str_curLine = linesInputDeck[iline].substr(0,ibackslash-1); // -1 since the backslash is not needed
        std::string str_nextLine = linesInputDeck[iline+1]; // the next line

        linesInputDeck[iline].clear(); // empty the current line
        linesInputDeck[iline+1] = str_curLine + str_nextLine; // merge into the next line

      } // iline
    } // merge
  } // iline



  // Step 3: push back non-trivial input-deck lines into tmp_vector
  int num_tmp_vector = 0;
  std::vector<std::string> tmp_vector;
  for (int iline = 0; iline < numLinesInputDeck; iline++) {
    if (linesInputDeck[iline].size() > 0) {

      num_tmp_vector++;
      tmp_vector.push_back(linesInputDeck[iline]);

    } // linesInputDeck[iline].size()
  } // iline
  //
  numLinesInputDeck = num_tmp_vector;
  linesInputDeck.clear();
  linesInputDeck = tmp_vector;
  tmp_vector.clear();
  //
  //str_counter << numLinesInputDeck;
  //MESSAGE_STDOUT(str_counter.str() + " lines exist the input deck.");
  //if (mpi::irank == 0)
  //  for (std::vector<std::string>::iterator it = linesInputDeck.begin(); it != linesInputDeck.end(); it++)
  //    std::cout << "After lines are merged: " << *it << "<EOL>" << std::endl;



  // Step 4: parse items line by line (see http://www.cplusplus.com/reference/cstring/strtok/)
  numEntriesInputDeck = numLinesInputDeck;
  entriesInputDeck = new t_EntryInputDeck[numEntriesInputDeck];
  //
  for (int iline = 0; iline < numLinesInputDeck; iline++) {

    char *cur_entry = (char *)cstr_to_constchar(linesInputDeck[iline]);
    char *token = NULL;
    token = strtok(cur_entry, delimiter); // strtok tokenizes a string using specified delimiters
    std::vector<std::string> tmp;

    while (token != NULL) {

      //if (mpi::irank == 0) std::cout << token << std::endl;
      tmp.push_back(std::string(token));
      token = strtok(NULL, delimiter);

    } // token
    // didn't allocate memory for cur_line and token; so no need to deallocate them

    // save parsed items
    entriesInputDeck[iline].name = tmp[FIRST];
    entriesInputDeck[iline].id = iline;
    entriesInputDeck[iline].number_items = tmp.size() - 1; // -1 since name is already stored above
    for (int i = SECOND; i < tmp.size(); i++)
      (entriesInputDeck[iline].body).push_back(tmp[i]);

    // check redundancy in name: if the same name shows up in the previous entries, remove the previous ones
    // exceptions (i.e. multiple instances of name are allowed)
    if (entriesInputDeck[iline].name == "PROBE")
      continue;
    for (int iline_prev = 0; iline_prev < iline; iline_prev++) {
      if (entriesInputDeck[iline_prev].name == entriesInputDeck[iline].name) {

        MESSAGE_STDOUT("A redundant input-deck item " + entriesInputDeck[iline].name + " is detected: take whatever shows up later.");

        entriesInputDeck[iline_prev].name.clear();
        entriesInputDeck[iline_prev].id = NONE;
        entriesInputDeck[iline_prev].number_items = 0;
        entriesInputDeck[iline_prev].body.clear();

      } // entriesInputDeck[iline_prev].name
    } // iline_prev
    tmp.clear();

  } // iline
  //if (mpi::irank == 0)
  //  for (int ientry = 0; ientry < numEntriesInputDeck; ientry++)
  //    if (entriesInputDeck[ientry].number_items > 0)
  //      std::cout << entriesInputDeck[ientry].name << std::endl;

  return;

} // parse_linesInputDeck



void clear_inputDeck(void) {

  numLinesInputDeck = 0;
  linesInputDeck.clear();

  numEntriesInputDeck = 0;
  delete[] entriesInputDeck;

  return;

} // clear_inputDeck



int count_inputDeck_name(std::string name) {

  // given name, count input-deck entries starting with this name
  int count = 0;
  for (int ientry = 0; ientry < numEntriesInputDeck; ientry++) {
    if (name == entriesInputDeck[ientry].name) {

      count++;

    } // name
  } // ientry

  return count;

} // count_inputDeck_name



int check_inputDeck_name(std::string name, int name_count) {

  // given name, find an <name_count>-th input-deck entry starting with this name and return its index in entriesInputDeck
  // useful if there exist multiple input-deck entries starting with this name
  int found = NONE;
  int counter = 0;
  for (int ientry = 0; ientry < numEntriesInputDeck; ientry++) {
    if (name == entriesInputDeck[ientry].name) {

      if (counter == name_count) {
        found = ientry;
        break;
      } // counter
      counter++;

    } // name
  } // ientry
  if (found == NONE)
    mpi::graceful_exit("The name " + name + " does not match any of the names in the input deck.");

  return found;

} // check_inputDeck_name



int check_inputDeck_keyword(std::string name, std::string keyword, int name_count) {

  // first get the index for name
  int index_name = check_inputDeck_name(name, name_count);

  // check if this name contains this keyword
  int found = NONE;
  for (int item = 0; item < entriesInputDeck[index_name].number_items; item++) {
    if (keyword == entriesInputDeck[index_name].body[item]) {

      found = item;
      break;

    } // keyword
  } // item
  if (found == NONE) {

    MESSAGE_STDOUT("The keyword " + keyword + " does not match any of the keywords in the list named as " + name + ".");
    MESSAGE_STDOUT("The list named as " + name + " has:");
    if (mpi::irank == 0)
      for (int item = 0; item < entriesInputDeck[index_name].number_items; item++)
        std::cout << entriesInputDeck[index_name].body[item] << std::endl;
    mpi::graceful_exit("The code terminates here.");

  } // found

  return found;

} // check_inputDeck_keyword



void get_userInput(std::string name, int &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name
  int index_name = check_inputDeck_name(name, name_count);

  std::string tmp = entriesInputDeck[index_name].body[0];

  // convert and store
  data = atoi(cstr_to_constchar(tmp));

  return;

} // get_userInput



void get_userInput(std::string name, int count, int *&data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name
  int index_name = check_inputDeck_name(name, name_count);

  // extract items right to this name "count" times
  assert( count > 0 );
  std::vector<std::string> tmp;
  for (int i = 0; i < count; i++)
    tmp.push_back(entriesInputDeck[index_name].body[i]);

  // convert and store
  for (int i = 0; i < count; i++)
    data[i] = atoi(cstr_to_constchar(tmp[i]));
  tmp.clear();

  return;

} // get_userInput



void get_userInput(std::string name, double &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name
  int index_name = check_inputDeck_name(name, name_count);

  std::string tmp = entriesInputDeck[index_name].body[0];

  // convert and store
  data = atof(cstr_to_constchar(tmp));

  return;

} // get_userInput



void get_userInput(std::string name, int count, double *&data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name
  int index_name = check_inputDeck_name(name, name_count);

  // extract items right to this name "count" times
  assert( count > 0 );
  std::vector<std::string> tmp;
  for (int i = 0; i < count; i++)
    tmp.push_back(entriesInputDeck[index_name].body[i]);

  // convert and store
  for (int i = 0; i < count; i++)
    data[i] = atof(cstr_to_constchar(tmp[i]));
  tmp.clear();

  return;

} // get_userInput



void get_userInput(std::string name, std::string &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name
  int index_name = check_inputDeck_name(name, name_count);

  data = entriesInputDeck[index_name].body[0];

  return;

} // get_userInput



void get_userInput(std::string name, int count, std::vector<std::string> &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name
  int index_name = check_inputDeck_name(name, name_count);

  data.clear();

  // extract items right to this name "count" times
  assert( count > 0 );
  for (int i = 0; i < count; i++)
    data.push_back(entriesInputDeck[index_name].body[i]);

  return;

} // get_userInput



void get_userInput(std::string name, std::string keyword, int &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name and keyword
  int index_name = check_inputDeck_name(name, name_count);
  int index_keyword = check_inputDeck_keyword(name, keyword, name_count);

  std::string tmp = entriesInputDeck[index_name].body[index_keyword+1]; // +1 since the keyword comes first

  // convert and store
  data = atoi(cstr_to_constchar(tmp));

  return;

} // get_userInput



void get_userInput(std::string name, std::string keyword, int count, int *&data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name and keyword
  int index_name = check_inputDeck_name(name, name_count);
  int index_keyword = check_inputDeck_keyword(name, keyword, name_count);

  // extract items right to this keyword "count" times
  assert( count > 0 );
  std::vector<std::string> tmp;
  for (int i = 0; i < count; i++)
    tmp.push_back(entriesInputDeck[index_name].body[index_keyword+1+i]); // +1 since the keyword comes first

  // convert and store
  for (int i = 0; i < count; i++)
    data[i] = atoi(cstr_to_constchar(tmp[i]));
  tmp.clear();

  return;

} // get_userInput



void get_userInput(std::string name, std::string keyword, double &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name and keyword
  int index_name = check_inputDeck_name(name, name_count);
  int index_keyword = check_inputDeck_keyword(name, keyword, name_count);

  std::string tmp = entriesInputDeck[index_name].body[index_keyword+1]; // +1 since the keyword comes first

  // convert and store
  data = atof(cstr_to_constchar(tmp));

  return;

} // get_userInput



void get_userInput(std::string name, std::string keyword, int count, double *&data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name and keyword
  int index_name = check_inputDeck_name(name, name_count);
  int index_keyword = check_inputDeck_keyword(name, keyword, name_count);

  // extract items right to this keyword "count" times
  assert( count > 0 );
  std::vector<std::string> tmp;
  for (int i = 0; i < count; i++)
    tmp.push_back(entriesInputDeck[index_name].body[index_keyword+1+i]); // +1 since the keyword comes first

  // convert and store
  for (int i = 0; i < count; i++)
    data[i] = atof(cstr_to_constchar(tmp[i]));
  tmp.clear();

  return;

} // get_userInput



void get_userInput(std::string name, std::string keyword, std::string &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name and keyword
  int index_name = check_inputDeck_name(name, name_count);
  int index_keyword = check_inputDeck_keyword(name, keyword, name_count);

  data = entriesInputDeck[index_name].body[index_keyword+1]; // +1 since the keyword comes first

  return;

} // get_userInput



void get_userInput(std::string name, std::string keyword, int count, std::vector<std::string> &data, int name_count) { // if name_count is not provided, name_count = 0 (see the header file)

  // get indices for name and keyword
  int index_name = check_inputDeck_name(name, name_count);
  int index_keyword = check_inputDeck_keyword(name, keyword, name_count);

  data.clear();

  // extract items right to this keyword "count" times
  assert( count > 0 );
  for (int i = 0; i < count; i++)
    data.push_back(entriesInputDeck[index_name].body[index_keyword+1+i]); // +1 since the keyword comes first

  return;

} // get_userInput

} // inputDeck
