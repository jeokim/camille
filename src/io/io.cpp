#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <locale>
#include <algorithm>
#include <cstddef>

#include "io.h"

namespace io {

// instantiate externs
ReportParallel report_parallel;
ReportGeometry report_geometry;
ReportUserInput report_userInput;
ReportInitialData report_initialData;
ReportSolution report_solution;
ReportSolution report_timing;

int counter_solution_files = 0;

t_BoundaryData *boundaryData;



void create_directory_basic() {

  // parallel
  report_parallel.path = path_report_parallel;
  mkdir(report_parallel.path);

  // geometry
  report_geometry.path = path_report_geometry;
  mkdir(report_geometry.path);

  // user input
  report_userInput.path = path_report_userInput;
  mkdir(report_userInput.path);

  // initial data
  report_initialData.path = path_report_initialData;
  mkdir(report_initialData.path);

  // solution
  report_solution.path = path_report_solution;
  mkdir(report_solution.path);

  // timing
  report_timing.path = path_report_timing;
  mkdir(report_timing.path);

  return;

} // create_directory_basic



void mkdir(std::string name_dir) {

  std::string command = "mkdir -p " + name_dir;

  if (mpi::irank == 0)
    system(cstr_to_constchar(command));

  return;

} // create_directory



void ReportParallel::initial_core() {

  std::ofstream ofs;

  for (int irank = 0; irank < mpi::nprocs; irank++) {
    if (mpi::irank == irank) { // write in the order of rank

      ofs.open(cstr_to_constchar(report_parallel.path + filename_parallel_core), std::ofstream::app);
      ofs << "Rank " << std::setw(8) << mpi::irank << std::setw(8) << "out of " << mpi::nprocs << " ready to run." << std::endl;
      ofs.close();

    } // mpi::irank
    mpi::wait_allothers();

  } // irank

  return;

} // ReportParallel::initial_core



void ReportParallel::initial_workLoad(Geometry::StructuredGrid *mygrid, int num_cells_total) {

  std::ofstream ofs;

  for (int irank = 0; irank < mpi::nprocs; irank++) {
    if (mpi::irank == irank) { // write in the order of rank

      ofs.open(cstr_to_constchar(report_parallel.path + filename_parallel_workLoad), std::ofstream::app);
      if (mpi::irank == 0)
        ofs << "variables = rank, block-ID, #-of-points-per-core, #-of-points-per-core(%), #-of-points-per-core-including-ghost-points" << std::endl;
      ofs << std::setw(8) << mpi::irank << "   "
          << std::setw(8) << mygrid->id_parent << "   "
          << std::setw(8) << mygrid->num_cells << "   "
          << std::scientific << static_cast<double>(mygrid->num_cells)/static_cast<double>(num_cells_total)*100.0 << "   "
          << std::setw(8) << mygrid->num_ocells << std::endl;
      ofs.close();

    } // mpi::irank
    mpi::wait_allothers();

  } // irank

  return;

} // ReportParallel::initial_workLoad



void ReportGeometry::initial_region(Geometry::Generic *region) {

  std::ofstream ofs;

  // region corresponds to an entire simulation domain; thus, every processor has the same information as to region
  if (mpi::irank == 0) {

    ofs.open(cstr_to_constchar(report_geometry.path + filename_geometry), std::ofstream::app);
    ofs << ">>> Region level" << std::endl;
    ofs << "Region has " << std::setw(4) << region->num_partitions << " blocks." << std::endl;
    ofs << std::endl;
    ofs.close();

  }

  mpi::wait_allothers();

  return;

} // ReportGeometry::initial



void ReportGeometry::initial_block_structured(Geometry::StructuredBlock *block) {

  std::ofstream ofs;

  // information on all blocks is shared by all cores
  if (mpi::irank == 0) {

    // file header
    ofs.open(cstr_to_constchar(report_geometry.path + filename_geometry), std::ofstream::app);
    ofs << ">>> Block level" << std::endl;
    ofs.close();

    // block ID and number of cells
    ofs.open(cstr_to_constchar(report_geometry.path + filename_geometry), std::ofstream::app);
    for (int iblock = 0; iblock < block[0].num_ourkinds; iblock++) {

      ofs << "ID: " << std::setw(4) << iblock << ", cells in XI: " << std::setw(6) << block[iblock].num_cells_dir[XI]   << 
                                                         ", ETA: " << std::setw(6) << block[iblock].num_cells_dir[ETA]  << 
                                                        ", ZETA: " << std::setw(6) << block[iblock].num_cells_dir[ZETA] << 
                                                       "; total: " << std::setw(6) << block[iblock].num_cells           << std::endl;
    } // iblock
    ofs << std::endl;
    ofs.close();

    //// block ID and number of nodes
    //ofs.open(cstr_to_constchar(report_geometry.path + filename_geometry), std::ofstream::app);
    //for (int iblock = 0; iblock < block[0].num_ourkinds; iblock++) {
    //
    //  ofs << "ID: " << std::setw(4) << iblock << ", nodes in XI: " << std::setw(6) << block[iblock].num_nodes_dir[XI]   << 
    //                                                     ", ETA: " << std::setw(6) << block[iblock].num_nodes_dir[ETA]  << 
    //                                                    ", ZETA: " << std::setw(6) << block[iblock].num_nodes_dir[ZETA] << 
    //                                                   "; total: " << std::setw(6) << block[iblock].num_nodes           << std::endl;
    //} // iblock
    //ofs << std::endl;
    //ofs.close();

    // block ID and number of cores and points
    ofs.open(cstr_to_constchar(report_geometry.path + filename_geometry), std::ofstream::app);
    for (int iblock = 0; iblock < block[0].num_ourkinds; iblock++) {

      ofs << "ID: " << std::setw(4) << iblock << ", cores (cells) in XI: " << std::setw(3) << block[iblock].num_cores_dir[XI]   << " ("<< std::setw(6) << block[iblock].num_cells_dir[XI]   << " )" <<
                                                                 ", ETA: " << std::setw(3) << block[iblock].num_cores_dir[ETA]  << " ("<< std::setw(6) << block[iblock].num_cells_dir[ETA]  << " )" <<
                                                                ", ZETA: " << std::setw(3) << block[iblock].num_cores_dir[ZETA] << " ("<< std::setw(6) << block[iblock].num_cells_dir[ZETA] << " )" << 
                                                               "; total: " << std::setw(3) << block[iblock].num_cores           << " ("<< std::setw(6) << block[iblock].num_cells           << " )" << std::endl;

    } // iblock
    ofs << std::endl;
    ofs.close();

  } // mpi::irank

  mpi::wait_allothers();

  return;

} // ReportGeometry::initial_block_structured



void ReportGeometry::initial_grid_structured(Geometry::StructuredGrid *grid) {

  std::ofstream ofs;

  // each grid belongs to different cores; thus, loop over entire cores

  // at the parent-level block
  for (int irank = 0; irank < mpi::nprocs; irank++) {

    if (mpi::irank == irank) { // write in the order of rank

      ofs.open(cstr_to_constchar(report_geometry.path + filename_geometry), std::ofstream::app);

      if (mpi::irank == 0)
        ofs << ">>> Grid level (parent)" << std::endl;

      ofs << "Global ID: "       << std::setw(3) << grid->id_global;
      ofs << ", local ID: "      << std::setw(3) << grid->id_local;
      ofs << ", parent block: "  << std::setw(3) << grid->id_parent;
      ofs << ", rank: "          << std::setw(4) << mpi::irank;

      ofs << ", is-ie: " << std::setw(4) << grid->is_in_parent[XI]    << ", " << std::setw(4) << grid->ie_in_parent[XI] <<
                    ", " << std::setw(4) << grid->is_in_parent[ETA]   << ", " << std::setw(4) << grid->ie_in_parent[ETA] <<
                    ", " << std::setw(4) << grid->is_in_parent[ZETA]  << ", " << std::setw(4) << grid->ie_in_parent[ZETA] << 
           ", iso-ieo: " << std::setw(4) << grid->iso_in_parent[XI]   << ", " << std::setw(4) << grid->ieo_in_parent[XI] <<
                    ", " << std::setw(4) << grid->iso_in_parent[ETA]  << ", " << std::setw(4) << grid->ieo_in_parent[ETA] <<
                    ", " << std::setw(4) << grid->iso_in_parent[ZETA] << ", " << std::setw(4) << grid->ieo_in_parent[ZETA] << std::endl;

      if (mpi::irank == mpi::nprocs - 1)
        ofs << std::endl;

      ofs.close();

    } // mpi::irank

    mpi::wait_allothers();

  } // irank

  mpi::wait_allothers();

  // at the local level
  for (int irank = 0; irank < mpi::nprocs; irank++) {

    if (mpi::irank == irank) { // write in the order of rank

      ofs.open(cstr_to_constchar(report_geometry.path + filename_geometry), std::ofstream::app);

      if (mpi::irank == 0)
        ofs << ">>> Grid level (local)" << std::endl;

      ofs << "Global ID: "       << std::setw(3) << grid->id_global;
      ofs << ", local ID: "      << std::setw(3) << grid->id_local;
      ofs << ", parent block: "  << std::setw(3) << grid->id_parent;
      ofs << ", rank: "          << std::setw(4) << mpi::irank;

      ofs << ", is-ie: " << std::setw(4) << grid->is[XI]    << ", " << std::setw(4) << grid->ie[XI] <<
                    ", " << std::setw(4) << grid->is[ETA]   << ", " << std::setw(4) << grid->ie[ETA] <<
                    ", " << std::setw(4) << grid->is[ZETA]  << ", " << std::setw(4) << grid->ie[ZETA] << 
           ", iso-ieo: " << std::setw(4) << grid->iso[XI]   << ", " << std::setw(4) << grid->ieo[XI] <<
                    ", " << std::setw(4) << grid->iso[ETA]  << ", " << std::setw(4) << grid->ieo[ETA] <<
                    ", " << std::setw(4) << grid->iso[ZETA] << ", " << std::setw(4) << grid->ieo[ZETA] << std::endl;

      if (mpi::irank == mpi::nprocs - 1)
        ofs << std::endl;

      ofs.close();

    } // mpi::irank

    mpi::wait_allothers();

  } // irank

  mpi::wait_allothers();

  return;

} // ReportGeometry::initial_grid_structured



void ReportUserInput::write_userInput(UserInput myinput) {

  std::ofstream ofs;

  if (mpi::irank == 0) {

    ofs.open(cstr_to_constchar(report_userInput.path + filename_userInput), std::ofstream::app);

    // geometry
    ofs << "Number of dimension: " << myinput.num_dim << std::endl;
    ofs << "Axisymmetry?: ";
    if (myinput.axisym == TRUE)
      ofs << "axisymmetric along the " << myinput.int_2xyz(myinput.axis_of_sym) << " axis" << std::endl;
    else
      ofs << "not axisymmetric" << std::endl;

    // physical model
    ofs << "Physical model: " << myinput.model_pde << std::endl;
    ofs << "Fluid model: " << myinput.model_fluid << std::endl;
    if (myinput.num_scalar > 0)
      ofs << "Number of passive scalars: " << myinput.num_scalar << std::endl;

    // simulation
    ofs << "Simulation: " << myinput.simulation << std::endl;

    // thermodynamics
    ofs << "Ratio of specific heats, gamma = c_p/c_v (if constant): " << myinput.gamma_specificheat << std::endl;
    ofs << "Specific heat at constant pressure, c_p (if constant): " << myinput.c_p << std::endl;

    // temporal discretization
    ofs << "Which quantity is fixed during time-stepping? (e.g. DT, CFL, ...): " << myinput.fix_dt_or_cfl << std::endl;
    if (myinput.fix_dt_or_cfl == "FIX_DT")
      ofs << "Time-step size is fixed as " << myinput.dt << std::endl;
    else if (myinput.fix_dt_or_cfl == "FIX_CFL")
      ofs << "CFL number is fixed as " << myinput.cfl << std::endl;
    ofs << "Time-advancement scheme: " << myinput.temporal_scheme << std::endl;
    ofs << "Number of total time steps: " << myinput.num_time_steps << std::endl;
    ofs << "Number of reporting time steps: " << myinput.report_freq << std::endl;

    // spatial discretization
    ofs << "Spatial discretization scheme: " << myinput.spatial_scheme << std::endl;
    if (myinput.spatial_scheme == "STANDARD_CENTRAL")
      ofs << "Order of accuracy: " << myinput.OA_spatial << std::endl;

    // low-pass filter
    if (myinput.do_filter == TRUE) {

      ofs << "Low-pass filter: " << myinput.filter_scheme << std::endl;
      if (myinput.filter_scheme == "STANDARD_CENTRAL") {

        ofs << "Filter order of accuracy: " << myinput.OA_filter << std::endl;

      } // myinput.filter_scheme
      else if (myinput.filter_scheme == "SFO11P") {

        ofs << "Filter strength: " << myinput.filter_strength << std::endl;

      } // myinput.filter_scheme
      ofs << "Filter-blend parameter: " << myinput.filter_blend << std::endl;

    } // myinput.do_filter
    else
      ofs << "No low-pass filter is applied; some sort of numerical stabilization would be necessary though." << std::endl;

    // buffer zone
    ofs << "If used, buffer-zone parameters are: " << "polynomial order: " << myinput.buffer_polynomial_order
                                                   << ", constant: " << myinput.buffer_constant << std::endl;

    // metric computation
    ofs << "Metric computation: " << myinput.scheme_metrics << std::endl;

    // overset
    if (myinput.do_overset == TRUE)
      ofs << "Overset-grid interpolation: " << myinput.overset_format << " format" << std::endl;
      if (myinput.overset_format == "OVERTURE") {

        ofs << "Overset interpolation type: " << myinput.overset_type_of_interpolation << std::endl;
        ofs << "Overset interpolation order of accuracy: " << myinput.overset_accuracy << std::endl;

      } // myinput.overset_format

    // domain decomposition
    if (myinput.how2decompDomain == "DECOMP_1D")
      ofs << "Domain is decomposed in the direction having the largest number of grid points." << std::endl;

    else if (myinput.how2decompDomain == "DECOMP_FROMFILE")

      ofs << "Domain is decomposed using information provided by " << myinput.file_decompDomain << std::endl;

    else {

      MESSAGE_STDOUT("Unknown domain decomposition!");
      assert(0);

    } // myinput.how2decompDomain

    // file
    ofs << "I/O file format: " << myinput.type_file << std::endl;
    if (myinput.type_file != "PLOT3D") {

      MESSAGE_STDOUT("PLOT3D is the only supported file format.");
      assert(0);

    } // myinput.type_file

    ofs << "You grid file: " << myinput.file_grid_in << std::endl;
    ofs << "You overset file: " << myinput.file_overset << std::endl;

    if (myinput.present_file_solution_in == TRUE)
      ofs << "Simulation restarts from " + myinput.file_solution_in << std::endl;
    else
      ofs << "Simulation starts from scratch" << std::endl;

    if (myinput.present_file_mean_in == TRUE)
      ofs << "Your base-state file: " << myinput.file_mean_in << std::endl;
    if (myinput.present_file_aux_in == TRUE)
      ofs << "A file containing auxiliary variables: " << myinput.file_aux_in << std::endl;

    ofs << "You boundary-condition file: " << myinput.file_boundary << std::endl;

    ofs << "Solution is written at every " << myinput.time_writing_solutions << " code time" << std::endl;

    if (myinput.harmonicWave_waveType != "NONE")
      ofs << "Time-harmonic wave is used in the simulation; amplitude: " << myinput.harmonicWave_amplitude <<
                                                             ", shape: " << myinput.harmonicWave_waveForm <<
                                             ", propagation direction: " << myinput.int_2xyz(myinput.harmonicWave_idir_propagation) <<
                                                        ", wavelength: " << myinput.harmonicWave_wavelength << std::endl;

    if (myinput.interp_fromWhichFormat != "NONE") {

      ofs << "Solution interpolation; format: " << myinput.interp_fromWhichFormat <<
                       ", solution dimension: " << myinput.num_dim_interpSource <<
                       ", number of variables to interpolate: " << myinput.num_vars_interpSource <<
                       ", variables to interpolate: " << myinput.vars_interpSource <<
                       ", number of zones in the source file: " << myinput.num_zones_interpSource << std::endl;
      ofs << "Scaling quantities for interpolated solutions; ";
      if (myinput.vars_interpSource == "PRHOU")
        ofs << "pressure: " << myinput.scale_p << "; density: " << myinput.scale_rho << "; velocity: " << myinput.scale_u << std::endl;
      else if (myinput.vars_interpSource == "PUS")
        ofs << "pressure: " << myinput.scale_p << "; velocity: " << myinput.scale_u << "; density: " << myinput.scale_s << std::endl;

      ofs << "Number of low-pass filters applied to interpolated solutions: " << myinput.num_filters_interpolation << std::endl;
      ofs << "Interpolated solutions are written in a format of " << myinput.interpolate_into_which_PLOT3D << std::endl;

    } // myinput.interp_fromWhichFormat

    ofs.close();

  } // mpi::irank

  mpi::wait_allothers();

  return;

}



void ReportSolution::minmax(int time_step, double time_sol, Geometry::StructuredGrid *mygrid, State *mystate) {

  std::ofstream ofs;

  int num_vars = mystate->num_vars_sol;

  double *sol_min = new double[num_vars];
  double *sol_max = new double[num_vars];
  double *sol_min_global = new double[num_vars];
  double *sol_max_global = new double[num_vars];
  for (int ivar = 0; ivar < num_vars; ivar++) {

    sol_min[ivar] =  abs(DUMMY_DOUBLE);
    sol_max[ivar] = -abs(DUMMY_DOUBLE);
    sol_min_global[ivar] = DUMMY_DOUBLE;
    sol_max_global[ivar] = DUMMY_DOUBLE;

  } // ivar

  // find the minimum and maximum values of solution variables at this core
  for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++)
    for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++)
      for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {

        int l0 = mygrid->idx1D(i, j, k);

        for (int ivar = 0; ivar < num_vars; ivar++) {

          sol_min[ivar] = std::min( sol_min[ivar], (mystate->sol[ivar])[l0] );
          sol_max[ivar] = std::max( sol_max[ivar], (mystate->sol[ivar])[l0] );

        } // ivar
      } // i
  mpi::wait_allothers();

  MPI_Allreduce(sol_min, sol_min_global, num_vars, MPI_DOUBLE, MPI_MIN, mpi::comm_region);
  MPI_Allreduce(sol_max, sol_max_global, num_vars, MPI_DOUBLE, MPI_MAX, mpi::comm_region);

  if (mpi::irank == 0) {

    ofs.open(cstr_to_constchar(report_solution.path + filename_solution), std::ofstream::app);
    ofs << "Time step: " << std::setw(10) << time_step
         << ", time: " << std::scientific << time_sol;
    for (int ivar = 0; ivar < num_vars; ivar++)
      ofs << ", " << mystate->name_vars[ivar] << ": " << std::scientific << sol_min_global[ivar] << " ~ "
                                                      << std::scientific << sol_max_global[ivar];
    ofs << std::endl;
    ofs.close();

  } // mpi::irank

  delete[] sol_min, sol_max, sol_min_global, sol_max_global;

  mpi::wait_allothers();

  return;

} // ReportSolution::minmax



void report_timeadvancing(int time_step, double time_sol, double dt, double cfl, UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate) {

  //if ( time_step%(myinput->report_freq) != 0)
  //  return;

  if (mpi::irank == 0)
    std::cout << "Time step: " << std::setw(10) << time_step
              << ", time: " << std::scientific << time_sol
              << ", dt: " << std::scientific
              << dt << ", CFL: " << std::scientific << cfl <<std::endl;

  report_solution.minmax(time_step, time_sol, mygrid, mystate);

  return;

} // report_timeadvancing



int skip_this_line_of_textFile(std::string line_cur) {

  int skip = FALSE; // by default, do not skip this line

  if (line_cur.size() == 0) // skip if an empty line
    skip = TRUE;

  else { // the string length is greater than zero
    std::string line_cur_copy = line_cur; // work on a copy of the original string
    line_cur_copy.erase(remove_if(line_cur_copy.begin(), line_cur_copy.end(), isspace), line_cur_copy.end()); // try removing all white spaces

    if (line_cur_copy.size() == 0) // so it was a bunch of whitespaces
      skip = TRUE;

    else { // one last check for commented lines

      // find out the first non-white-space character since the comment symbol is not always at the first column
      int i = 0;
      int keep_going = TRUE;
      while( keep_going ) {
        char c = line_cur.at(i);
        if (isspace(c)) // use isspace instead of comparing c to ' '; whitespace can be ' ', '\t', '\n', ...
          i++;

        else {
          keep_going = FALSE; // if non-white-space, escape
          break;
        } // ~isspace(c)

        if (i == line_cur.size()) keep_going = FALSE; // if the end of the string is reached, exit anyway
      } // keep_going
      assert( i < line_cur.size() );
      //
      if (line_cur.at(i) == '#') // skip if a commented line
        skip = TRUE;

    }
  } // line_cur.size()

  return skip;

} // skip_this_line_of_textFile



void read_inputDeck(int argc, char * argv[]) {

  std::ifstream ifs;
  std::string line_cur;
  std::stringstream str_counter;

  MESSAGE_STDOUT("Reading the input deck " + filename_inputDeck);

  int num_lines_inputDeck = 0;

  // count how many lines should be read
  ifs.open(cstr_to_constchar(filename_inputDeck), std::ifstream::in);
  if (!ifs.is_open())
    mpi::graceful_exit("The input deck " + filename_inputDeck + " does not exist.");

  // count the number of lines which are neither comments nor empty
  std::getline(ifs, line_cur);
  while (!ifs.eof()) { // read until end-of-file is reached
    if (skip_this_line_of_inputDeck(line_cur) == FALSE) { // count lines only if they are not skipped

      num_lines_inputDeck++;

    } // skip_this_line_of_inputDeck(line_cur)

    // read the next line whether the current line is skipped or not
    std::getline(ifs, line_cur);

  } // !ifs.eof()
  ifs.close();
  //
  assert ( num_lines_inputDeck > 0 );
  //str_counter << num_lines_inputDeck;
  //MESSAGE_STDOUT(str_counter.str() + " lines are read from the input deck.");

  // store input-deck lines
  inputDeck::numLinesInputDeck = num_lines_inputDeck;
  //
  num_lines_inputDeck = 0;
  ifs.open(cstr_to_constchar(filename_inputDeck), std::ifstream::in);
  std::getline(ifs, line_cur);
  while (!ifs.eof()) { // read until end-of-file is reached
    if (skip_this_line_of_inputDeck(line_cur) == FALSE) { // read lines only if they are not skipped

      inputDeck::linesInputDeck.push_back(line_cur);
      //if (mpi::irank == 0) std::cout << "In the input deck: " << inputDeck::linesInputDeck[num_lines_inputDeck] << std::endl;
      num_lines_inputDeck++;

    } // skip_this_line_of_inputDeck(line_cur)

    // read the next line whether the current line is skipped or not
    std::getline(ifs, line_cur);

  } // !ifs.eof()
  ifs.close();
  assert(num_lines_inputDeck == inputDeck::numLinesInputDeck);

  return;

} // read_inputDeck



int skip_this_line_of_inputDeck(std::string line_cur) {

  int skip = skip_this_line_of_textFile(line_cur);

  return skip;

} // skip_this_line_of_inputDeck



void read_grid_size(UserInput *myinput, int &num_blocks_in_file, int *&num_cells_in) {

  std::string filename = myinput->file_grid_in;

  if (myinput->present_file_grid_in) {

    if (myinput->type_file == "PLOT3D")
      plot3d::read_grid_header(filename, num_blocks_in_file, num_cells_in);

    else
      mpi::graceful_exit("Unknown type of grid files.");

  } // myinput->present_file_grid_in
  else {

    // GRID_HARDWIRED
    num_blocks_in_file = 2; // number of blocks within region
    ALLOCATE1D_INT_2ARG(num_cells_in, num_blocks_in_file, DIM_MAX);
    num_cells_in[0*DIM_MAX+XI] = 80; num_cells_in[0*DIM_MAX+ETA] = 24; num_cells_in[0*DIM_MAX+ZETA] = 1;
    num_cells_in[1*DIM_MAX+XI] = 40; num_cells_in[1*DIM_MAX+ETA] = 12; num_cells_in[1*DIM_MAX+ZETA] = 1;

  } // myinput->present_file_grid_in

  return;

} // read_grid_size



void read_grid(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  if (myinput->present_file_grid_in) {

    if (myinput->type_file == "PLOT3D")
      plot3d::read_grid_serialIO(myinput, mygrid, block);

    else
      mpi::graceful_exit("Unknown type of grid files.");

  } // myinput->present_file_grid_in
  else
    mygrid->hardwire_gridPoint(myinput, &block[mygrid->id_parent]);

  return;

} // read_grid



void write_grid(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  if (myinput->type_file == "PLOT3D")
    plot3d::write_grid_serialIO(myinput, mygrid, block);

  else
    mpi::graceful_exit("Unknown type of grid files.");

  return;

} // write_grid



void read_solution(std::string filename, UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate) {

  int num_vars = myinput->num_vars_sol;

  if (myinput->type_file == "PLOT3D")
    if (num_vars == DIM_MAX+2) // if 5 solution variables, go with the PLOT3D solution format
      plot3d::read_solution_serialIO(filename, myinput, mygrid, block, mystate);

    else // if more than 5 solution variables, need to use the PLOT3D function format
      plot3d::read_function_serialIO(filename, myinput, mygrid, block, num_vars, mystate->sol);

  else
    mpi::graceful_exit("Unknown type of solution files.");

  return;

} // read_solution



void read_meanState(std::string filename, UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, int num_vars, double **var) {

  if (myinput->type_file == "PLOT3D")
    plot3d::read_function_serialIO(filename, myinput, mygrid, block, num_vars, var);

  else
    mpi::graceful_exit("Unknown type of mean-state files.");

  return;

} // read_meanState



void read_auxvar(std::string filename, UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, int num_vars, double **var) {

  if (myinput->type_file == "PLOT3D")
    plot3d::read_function_serialIO(filename, myinput, mygrid, block, num_vars, var);

  else
    mpi::graceful_exit("Unknown type of auxiliary-variable files.");

  return;

} // read_auxvar



void write_solution(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate) {

  std::string filename;
  std::stringstream str_counter;
  int num_vars = myinput->num_vars_sol;

  if (myinput->type_file == "PLOT3D") {

    // set the file name
    str_counter << std::setw(6) << std::setfill('0') << counter_solution_files;
    filename = myinput->file_solution + "." + str_counter.str() + ".q";

    if (num_vars == DIM_MAX+2) { // if 5 solution variables, go with the PLOT3D solution format

      plot3d::write_solution_serialIO(myinput, mygrid, block, mystate, filename);
      plot3d::write_solution_namefile(myinput->file_varname_solution, mystate);

    } // num_vars
    else { // if more than 5 solution variables, need to use the PLOT3D function format

      double **func_2write = new double *[num_vars];

      for (int ivar = 0; ivar < num_vars; ivar++)
        func_2write[ivar] = new double[mygrid->num_ocells];

      for (int ivar = 0; ivar < num_vars; ivar++)
        for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
          (func_2write[ivar])[l0] = (mystate->sol[ivar])[l0];

      plot3d::write_function_serialIO(myinput, mygrid, block, num_vars, func_2write, filename);
      plot3d::write_function_namefile(myinput->file_varname_solution, num_vars, mystate->name_vars);

      DEALLOCATE_2DPTR(func_2write, num_vars);

    }

  } // myinput->type_file
  else
    mpi::graceful_exit("Unknown type of data files.");

  return;

} // write_solution



void write_solution_mean(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate) {

  int num_vars = mystate->num_vars_mean;
  double **func_2write = new double *[num_vars];
  std::string filename;

  if (num_vars == 0)
    return;

  if (myinput->type_file == "PLOT3D") {

    for (int ivar = 0; ivar < num_vars; ivar++)
      func_2write[ivar] = new double[mygrid->num_ocells];

    for (int ivar = 0; ivar < num_vars; ivar++)
      for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
        (func_2write[ivar])[l0] = (mystate->sol_mean[ivar])[l0];

    filename = "mean_" + myinput->file_solution + ".q";
    plot3d::write_function_serialIO(myinput, mygrid, block, num_vars, func_2write, filename);
    filename = "mean_" + myinput->file_varname_solution;
    plot3d::write_function_namefile(filename, num_vars, mystate->name_vars_mean);

    DEALLOCATE_2DPTR(func_2write, num_vars);

  } // myinput->type_file
  else
    mpi::graceful_exit("Unknown type of data files.");

  return;

} // write_solution_mean



void write_solution_aux(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate) {

  int num_vars = mystate->num_vars_aux;
  double **func_2write = new double *[num_vars];
  std::string filename;

  if (num_vars == 0)
    return;

  // solution counter is shared with solution files
  std::stringstream str_counter;
  str_counter << std::setw(6) << std::setfill('0') << counter_solution_files;

  if (myinput->type_file == "PLOT3D") {

    for (int ivar = 0; ivar < num_vars; ivar++)
      func_2write[ivar] = new double[mygrid->num_ocells];

    for (int ivar = 0; ivar < num_vars; ivar++)
      for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
        (func_2write[ivar])[l0] = (mystate->sol_aux[ivar])[l0];

    filename = "aux_" + myinput->file_solution + "." + str_counter.str() + ".q";
    plot3d::write_function_serialIO(myinput, mygrid, block, num_vars, func_2write, filename);
    filename = "aux_" + myinput->file_varname_solution;
    plot3d::write_function_namefile(filename, num_vars, mystate->name_vars_aux);

    DEALLOCATE_2DPTR(func_2write, num_vars);

  } // myinput->type_file
  else
    mpi::graceful_exit("Unknown type of data files.");

  return;

} // write_solution_aux



void write_solution_on_the_fly(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate, int num_time_steps) {

  if ( static_cast<int>( (mystate->time_sol - mystate->time_sol_lastrun) / myinput->time_writing_solutions ) > counter_solution_files ) {

    counter_solution_files++;
    write_solution(myinput, mygrid, block, mystate);

    // if there are any auxiliary variables, write them too
//    if (myinput->num_vars_aux > 0)
//      write_solution_aux(myinput, mygrid, block, mystate);

  } // static_cast<int>
  else if ( mystate->time_step == num_time_steps ) {

    counter_solution_files++;
    write_solution(myinput, mygrid, block, mystate);

    // if there are any auxiliary variables, write them too
//    if (myinput->num_vars_aux > 0)
//      write_solution_aux(myinput, mygrid, block, mystate);

  } // mystate->time_step

  mpi::wait_allothers();

  return;

}



void write_metrics(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  int num_vars = DIM_MAX * DIM_MAX + 2; // 3 by 3 metric tensor, Jacobian, and inverse Jacobian
  double **func_2write = new double *[num_vars];
  std::string filename;

  if (myinput->type_file == "PLOT3D") {

    for (int ivar = 0; ivar < num_vars; ivar++)
      func_2write[ivar] = new double[mygrid->num_ocells];

    for (int l0 = 0; l0 < mygrid->num_ocells; l0++) {

      int index_var = 0;
      for (int irow = XI; irow < DIM_MAX; irow++)
        for (int jcol = XDIR; jcol < DIM_MAX; jcol++) {
          func_2write[index_var][l0] = mygrid->cell[l0].metrics[irow][jcol];
          index_var++;
        } // jcol

      func_2write[DIM_MAX * DIM_MAX][l0] = mygrid->cell[l0].Jac;
      func_2write[DIM_MAX * DIM_MAX + 1][l0] = mygrid->cell[l0].invJac;

    } // l0

    filename = myinput->file_metrics;
    plot3d::write_function_serialIO(myinput, mygrid, block, num_vars, func_2write, filename);
    plot3d::write_function_namefile(myinput->file_varname_metrics, num_vars, plot3d::name_metrics());

    DEALLOCATE_2DPTR(func_2write, num_vars);

  } // myinput->type_file
  else
    mpi::graceful_exit("Unknown type of data files.");

  return;

} // write_metrics



void read_overset(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  if (myinput->overset_format == "OVERTURE_NOTSCALING")
    overture::read_overset_NOTSCALING(myinput, mygrid, block);

  else if (myinput->overset_format == "OVERTURE")
    overture::read_overset(myinput, mygrid, block);

  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  mpi::wait_allothers();

  return;

} // read_overset



void cleanup_overset_temporarydata(UserInput *myinput) {

  if (myinput->overset_format == "OVERTURE_NOTSCALING" || myinput->overset_format == "OVERTURE")
    overture::cleanup_overset_temporarydata(myinput);

  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  mpi::wait_allothers();

  return;

} // cleanup_overset_temporarydata



void read_bc(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  std::ifstream ifs;
  std::string line_cur;

  int num_boundaryCondition_nonperiodic = 0;
  int num_bufferZones = 0;

  t_BoundaryData tmpBoundaryData;

  int iso_in_block[DIM_MAX], ieo_in_block[DIM_MAX];
  int iso_in_grid[DIM_MAX], ieo_in_grid[DIM_MAX];

  // read the boundary condition file line by line to know how many boundaries there are
  ifs.open(cstr_to_constchar(myinput->file_boundary), std::ifstream::in);
  if (!ifs.is_open())
    mpi::graceful_exit("The file for your boundary treatment does not exist.");

  // read the boundary condition file
  std::getline(ifs, line_cur);
  while (!ifs.eof()) { // read until end-of-file is reached
    if (skip_this_line_of_bc_file(line_cur) == FALSE) { // take a look at each line only if we don't skip it

      extract_bcData_out_of_(line_cur, mygrid, block, tmpBoundaryData);

//      if (mpi::irank == 0) {
//        std::cout << line_cur << std::endl;
//        std::cout << "iblock: " << std::setw(4) << tmpBoundaryData.iblock
//                  << ", itype: " << std::setw(4) << tmpBoundaryData.itype
//                  << ", idir: " << std::setw(4) << tmpBoundaryData.idir
//                  << ", iend: " << std::setw(4) << tmpBoundaryData.iend
//                  << ", is~ie in XI: " << std::setw(6) << tmpBoundaryData.is[XI] << "~" << std::setw(6) << tmpBoundaryData.ie[XI]
//                  << ", is~ie in ETA: " << std::setw(6) << tmpBoundaryData.is[ETA] << "~" << std::setw(6) << tmpBoundaryData.ie[ETA]
//                  << ", is~ie in ZETA: " << std::setw(6) << tmpBoundaryData.is[ZETA] << "~" << std::setw(6) << tmpBoundaryData.ie[ZETA] << std::endl;
//        std::cout << std::endl;
//      } // mpi::irank

      int whether_this_is_my_boundary = check_if_mygrid_has_this_boundary(mygrid, tmpBoundaryData);
      switch ( whether_this_is_my_boundary ) {
      case FALSE: // this boundary has nothing to do with my grid; do nothing

        break;

      case BOUNDARY_BC:

        num_boundaryCondition_nonperiodic++;

        break;

      case BOUNDARY_BUFFERZONE:

        num_bufferZones++;

        break;

      default:

        mpi::graceful_exit("The code supports either boundary condition or buffer zone.");

        break;

      } // whether_this_is_my_boundary
    } // skip_this_line_of_bc_file(line_cur)

    // read the next line whether the current line is skipped or not
    std::getline(ifs, line_cur);

  } // !ifs.eof()

  // initialize
  mygrid->initialize_boundaryCondition(num_boundaryCondition_nonperiodic);
  mygrid->initialize_bufferZone(num_bufferZones);
  //
//  std::cout << "Rank: " << std::setw(4) << mpi::irank << ", block: " << std::setw(4) << mygrid->id_parent
//            << ", num_boundaryCondition_nonperiodic: " << std::setw(4) << mygrid->num_boundaryCondition_nonperiodic
//            << ", num_bufferZones: " << std::setw(4) << mygrid->num_bufferZones << std::endl;

  // read the boundary condition file once again to set up specific boundary data
  // first go back to the begning
  ifs.clear();
  ifs.seekg(0, ifs.beg);

  // count boundaries once again
  num_boundaryCondition_nonperiodic = 0;
  num_bufferZones = 0;

  // read the boundary condition file once again
  std::getline(ifs, line_cur);
  while (!ifs.eof()) { // read until end-of-file is reached
    if (skip_this_line_of_bc_file(line_cur) == FALSE) { // take a look at each line only if we don't skip it

      extract_bcData_out_of_(line_cur, mygrid, block, tmpBoundaryData);

      int whether_this_is_my_boundary = check_if_mygrid_has_this_boundary(mygrid, tmpBoundaryData);
      switch ( whether_this_is_my_boundary ) {
      case FALSE: // this boundary has nothing to do with my grid; do nothing

        break;

      case BOUNDARY_BC: case BOUNDARY_BUFFERZONE:

        for (int idir = XI; idir < DIM_MAX; idir++) {

          // take a portion of the boundary data belonging to my grid
          iso_in_block[idir] = std::max( mygrid->iso_in_parent[idir], tmpBoundaryData.is[idir] );
          ieo_in_block[idir] = std::min( mygrid->ieo_in_parent[idir], tmpBoundaryData.ie[idir] );

          iso_in_grid[idir] = iso_in_block[idir] - mygrid->iso_in_parent[idir];
          ieo_in_grid[idir] = ieo_in_block[idir] - mygrid->iso_in_parent[idir];

        } // idir

        if (whether_this_is_my_boundary == BOUNDARY_BC) {

          (mygrid->boundaryCondition[num_boundaryCondition_nonperiodic]).initialize(mygrid->num_dim, iso_in_grid, ieo_in_grid, tmpBoundaryData.idir, tmpBoundaryData.iend, tmpBoundaryData.itype);
          num_boundaryCondition_nonperiodic++;

        } // whether_this_is_my_boundary
        else if (whether_this_is_my_boundary == BOUNDARY_BUFFERZONE) {

          (mygrid->bufferZone[num_bufferZones]).initialize(mygrid->num_dim, iso_in_grid, ieo_in_grid, tmpBoundaryData.idir, tmpBoundaryData.iend, tmpBoundaryData.itype);
          (mygrid->bufferZone[num_bufferZones]).set_parameters_4bufferZone(myinput, tmpBoundaryData.is, tmpBoundaryData.ie);
          num_bufferZones++;

        } // whether_this_is_my_boundary

        break;

      default:

        mpi::graceful_exit("The code supports either boundary condition or buffer zone.");

        break;

      } // whether_this_is_my_boundary
    } // skip_this_line_of_bc_file(line_cur)

    // read the next line whether the current line is skipped or not
    std::getline(ifs, line_cur);

  } // !ifs.eof()
  assert ( num_boundaryCondition_nonperiodic == mygrid->num_boundaryCondition_nonperiodic );
  assert ( num_bufferZones == mygrid->num_bufferZones );

  ifs.close();



  for (int ibc = FIRST; ibc < num_boundaryCondition_nonperiodic; ibc++) {


  } // ibc
mpi::wait_allothers();
mpi::graceful_exit("Bye now!");

  return;

} // read_bc



int skip_this_line_of_bc_file(std::string line_cur) {

  int skip = skip_this_line_of_textFile(line_cur);

  return skip;

} // skip_this_line_of_bc_file



void extract_bcData_out_of_(std::string line_cur, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, t_BoundaryData &tmpBoundaryData) {

  const int num_entries_per_line = 9; // number of entries per line of the file

  std::string entry; // an extracted entry out of line_cur
  std::string sub_entry;
  int ientry = 0; // counts entries on this line
  int size_entry;
  int entry_in_integer; // entry converted to an integer value

  std::size_t loc_start = FIRST; // start to search from the first column
  std::size_t loc_end;

  // reset the temporary boundary data
  tmpBoundaryData.iblock = NONE;
  tmpBoundaryData.itype = NONE;
  tmpBoundaryData.idir = NONE;
  tmpBoundaryData.iend = NONE;
  for (int idir = XI; idir < DIM_MAX; idir++) {

    tmpBoundaryData.is[idir] = 0;
    tmpBoundaryData.ie[idir] = 0;

  } // idir

  int cur_block;
  int idir_cur;
  int num_cells_in_block_dir[DIM_MAX];

  while ( loc_start < line_cur.size() && ientry < num_entries_per_line ) {

    // extract one entry at a time
    loc_start = line_cur.find_first_not_of(' ', loc_start);
    loc_end = line_cur.find_first_of(' ', loc_start) - 1;
    entry = line_cur.substr(loc_start, loc_end - loc_start + 1);
    size_entry = entry.size();
    loc_start = loc_end + 1;

    // always check more than necessary
    if (entry.at(FIRST) == '#')
      mpi::graceful_exit("Some entry in the boundary condition file has a comment sign.");

    // parse each entry
    switch ( ientry ) {
    case FIRST: // block ID

      tmpBoundaryData.iblock = atoi(cstr_to_constchar(entry));

      cur_block = tmpBoundaryData.iblock;
      for (int idir = XI; idir < DIM_MAX; idir++)
        num_cells_in_block_dir[idir] = block[cur_block].num_cells_dir[idir];

      break;

    case SECOND: // type of boundary

      tmpBoundaryData.itype = parse_boundary_type(entry);

      break;

    case THIRD: // direction and side

      // which side is this boundary (LEFT/RIGHT)
      if (entry.at(FIRST) == '+')
        tmpBoundaryData.iend = LEFT;

      else if (entry.at(FIRST) == '-')
        tmpBoundaryData.iend = RIGHT;

      else
        mpi::graceful_exit("Unknown side for a boundary in the boundary condition file.");

      // get the direction
      if (size_entry == 3) { // since +XI or -XI has 3 letters

        sub_entry = entry.substr(SECOND, size_entry - 1); // omit the sign
        str_tolower(sub_entry); // convert to lower cases to reduce the number of if's below
        if (sub_entry == "xi")
          tmpBoundaryData.idir = XI;

      } // size_entry
      else if (size_entry == 4) { // since +ETA or -ETA has 4 letters

        sub_entry = entry.substr(SECOND, size_entry - 1); // omit the sign
        str_tolower(sub_entry); // convert to lower cases to reduce the number of if's below
        if (sub_entry == "eta")
          tmpBoundaryData.idir = ETA;

      } // size_entry
      else if (size_entry == 5) { // since +ZETA or -ZETA has 5 letters

        sub_entry = entry.substr(SECOND, size_entry - 1); // omit the sign
        str_tolower(sub_entry); // convert to lower cases to reduce the number of if's below
        if (sub_entry == "zeta")
          tmpBoundaryData.idir = ZETA;

      } // size_entry
      else
        mpi::graceful_exit("Unknown direction for a boundary in the boundary condition file.");

      break;

    case FOURTH: case SIXTH: case EIGHTH: // is1 or is2 or is3

      idir_cur = (ientry - FOURTH) / 2;

      sub_entry = entry;
      if (entry.at(FIRST) == '+' || entry.at(FIRST) == '-')
        sub_entry = entry.substr(SECOND, size_entry - 1); // extract only a number

      entry_in_integer = atoi(cstr_to_constchar(sub_entry));

      // correct if a negative index
      if (entry.at(FIRST) == '-')
        entry_in_integer = (num_cells_in_block_dir[idir_cur] - 1) - entry_in_integer;

      tmpBoundaryData.is[idir_cur] = entry_in_integer;

      break;

    case FIFTH: case SEVENTH: case NINTH: // ie1 or ie2 or ie3

      idir_cur = (ientry - FIFTH) / 2;

      sub_entry = entry;
      if (entry.at(FIRST) == '+' || entry.at(FIRST) == '-')
        sub_entry = entry.substr(SECOND, size_entry - 1); // extract only a number

      entry_in_integer = atoi(cstr_to_constchar(sub_entry));

      // correct if a negative index
      if (entry.at(FIRST) == '-')
        entry_in_integer = (num_cells_in_block_dir[idir_cur] - 1) - entry_in_integer;

      tmpBoundaryData.ie[idir_cur] = entry_in_integer;

      break;

    default:

      mpi::graceful_exit("Some entry in the boundary condition file is incorrect; check the file.");

      break;

    } // ientry
    ientry++;

  } // keep_searching_this_line
  if (ientry != num_entries_per_line)
    mpi::graceful_exit("Number of entries in the boundary condition file is incorrect.");

  return;

} // extract_bcData_out_of_



int parse_boundary_type(std::string str_in) {

  int itype;

  str_tolower(str_in); // convert to a lower case so that a comparison can be made with fewer number of if's

  if (str_in  == "dirichlet")
    itype = BC_DIRICHLET;

  else if (str_in  == "dirichlet_allzero")
    itype = BC_DIRICHLET_ALLZERO;

  else if (str_in  == "dirichlet_harmonicwave")
    itype = BC_DIRICHLET_HARMONICWAVE;

  else if (str_in  == "dirichlet_file")
    itype = BC_DIRICHLET_FILE;

  else if (str_in  == "neumann")
    itype = BC_NEUMANN;

  else if (str_in  == "wall_slip_kinematic")
    itype = BC_WALL_SLIP_KINEMATIC;

  else if (str_in  == "wall_slip_kinematic_x")
    itype = BC_WALL_SLIP_KINEMATIC_X;

  else if (str_in  == "wall_slip_kinematic_y")
    itype = BC_WALL_SLIP_KINEMATIC_Y;

  else if (str_in  == "wall_slip_kinematic_z")
    itype = BC_WALL_SLIP_KINEMATIC_Z;

  else if (str_in  == "wall_slip_kinematic_r")
    itype = BC_WALL_SLIP_KINEMATIC_Y;

  else if (str_in == "centerline_cart_norm2x")
    itype = BC_CENTERLINE_CART_NORM2X;

  else if (str_in == "centerline_cart_norm2y")
    itype = BC_CENTERLINE_CART_NORM2Y;

  else if (str_in == "centerline_cart_norm2z")
    itype = BC_CENTERLINE_CART_NORM2Z;

  else if (str_in == "centerline_axisym")
    itype = BC_CENTERLINE_AXISYM;

  else if (str_in  == "sponge_freund_ambient")
    itype = SPONGE_FREUND_AMBIENT;

  else if (str_in  == "sponge_freund_dirichlet")
    itype = SPONGE_FREUND_DIRICHLET;

  else if (str_in  == "sponge_freund_harmonicwave")
    itype = SPONGE_FREUND_HARMONICWAVE;

  else
    mpi::graceful_exit("Unknown type of boundary treatment in the boundary condition file.");

  return itype;

} // parse_boundary_type



int check_if_mygrid_has_this_boundary(Geometry::StructuredGrid *mygrid, t_BoundaryData tmpBoundaryData) {

  // check if this boundary belongs, at least, to my block
  if (mygrid->id_parent != tmpBoundaryData.iblock)
    return FALSE;

  // compare the index range between the boundary data and my grid
  int counter = 0;
  for (int idir = XI; idir < DIM_MAX; idir++) {

    if (mygrid->is_in_parent[idir] <= tmpBoundaryData.ie[idir] &&
        mygrid->ie_in_parent[idir] >= tmpBoundaryData.is[idir])
      counter++;

  } // idir
  if (counter != DIM_MAX)
    return FALSE;

  // at this point, my grid has at least a portion of this boundary
  switch ( tmpBoundaryData.itype ) {
  case BC_DIRICHLET: case BC_DIRICHLET_ALLZERO: case BC_DIRICHLET_HARMONICWAVE: case BC_DIRICHLET_FILE: // if a boundary condition
  case BC_NEUMANN:
  case BC_WALL_SLIP_KINEMATIC:
  case BC_WALL_SLIP_KINEMATIC_X: case BC_WALL_SLIP_KINEMATIC_Y: case BC_WALL_SLIP_KINEMATIC_Z:
  case BC_CENTERLINE_CART_NORM2X: case BC_CENTERLINE_CART_NORM2Y: case BC_CENTERLINE_CART_NORM2Z: case BC_CENTERLINE_AXISYM:

    return BOUNDARY_BC;

    break;

  case SPONGE_FREUND_AMBIENT: case SPONGE_FREUND_DIRICHLET: case SPONGE_FREUND_HARMONICWAVE: // if a buffer zone

    return BOUNDARY_BUFFERZONE;

    break;

  default:
    mpi::graceful_exit("Unknown type of boundary treatment while checking if this boundary belongs to my grid.");

  } // tmpBoundaryData.itype

  return FALSE;

} // check_if_mygrid_has_this_boundary



void read_file_decompDomain(std::string file_decompDomain, int iblock, int num_cores_dir[DIM_MAX]) {

  std::ifstream ifs;
  std::string line_cur;
  int found = FALSE;

  ifs.open(cstr_to_constchar(file_decompDomain), std::ifstream::in);
  if (!ifs.is_open())
    mpi::graceful_exit("The file for your domain decomposition does not exist.");

  // read the domain decomposition file
  std::getline(ifs, line_cur);
  while (!ifs.eof() && found == FALSE) { // read until end-of-file is reached, or 
                                        // until the decomposition information for my block is found
    if (skip_this_line_of_bc_file(line_cur) == FALSE) { // take a look at each line only if we don't skip it

      const int num_entries_per_line = 4; // number of entries per line of the file
                                          // column 0: block ID starting from 0
                                          // column 1: number of cores in the \xi direction in this block
                                          // column 2: number of cores in the \eta direction in this block
                                          // column 3: number of cores in the \zeta direction in this block
      std::string entry; // an extracted entry out of line_cur
      int ientry = 0; // counts entries on this line
      int size_entry;
      std::size_t loc_start = FIRST; // start to search from the first column
      std::size_t loc_end;

      // extract the first entry for the block ID
      loc_start = line_cur.find_first_not_of(' ', loc_start);
      loc_end = line_cur.find_first_of(' ', loc_start) - 1;
      entry = line_cur.substr(loc_start, loc_end - loc_start + 1);
      size_entry = entry.size();
      loc_start = loc_end + 1;

      if (atoi(cstr_to_constchar(entry)) == iblock) {

        ientry = 1;
        int idir = 0;

        while ( loc_start < line_cur.size() && ientry < num_entries_per_line ) {

          // extract one entry at a time
          loc_start = line_cur.find_first_not_of(' ', loc_start);
          loc_end = line_cur.find_first_of(' ', loc_start) - 1;
          entry = line_cur.substr(loc_start, loc_end - loc_start + 1);
          size_entry = entry.size();
          loc_start = loc_end + 1;

          num_cores_dir[idir++] = atoi(cstr_to_constchar(entry));

          ientry++;

        } // loc_start
        found = TRUE;
//        std::cout << "Rank: " << mpi::irank << ", block ID: " << iblock
//                                            << ", core #: "   << num_cores_dir[XI] << " (in XI), "
//                                                              << num_cores_dir[ETA] << " (in ETA), "
//                                                              << num_cores_dir[ZETA] << " (in ZETA)" << std::endl;

      } // atoi(cstr_to_constchar(entry))
    } // skip_this_line_of_bc_file(line_cur)

    // read the next line whether the current line is skipped or not
    std::getline(ifs, line_cur);

  } // !ifs.eof()
  ifs.close();

  if (found == FALSE)
    mpi::graceful_exit("Some block did not find the decomposition information in the given file.");

  return;

} // read_file_decompDomain

} // io
