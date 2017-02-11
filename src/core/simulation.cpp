#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>

#include "simulation.h"

namespace simulation {

// input
UserInput myinput;

// geometry
int num_region = 1;
int num_blocks_global = 0;
int *num_cells_in;
//
Geometry::Generic *region;
Geometry::StructuredBlock *block;
Geometry::StructuredGrid *grid;

// state
State state;



void initialize(int argc, char * argv[]) {

  // start MPI
  MPI_Init(&argc, &argv);
  mpi::init_parallel_global();



  // reporting file I/O
  io::create_directory_basic();
  io::report_parallel.initial_core();
  todos::add("Make a directory keeping all benchmark problems (grid, input file, overset data, boundary condition file, ...)");



  // user inputs
  io::read_inputDeck(argc, argv);
  inputDeck::parse_linesInputDeck();
  myinput.set(argc, argv);
  io::report_userInput.write_userInput(myinput);
  inputDeck::clear_inputDeck();
  mpi::wait_allothers("User input is initialized.");



  // get minimal information on grid (total numbers of blocks, total numbers of cells, number of cells per each block, and so on)
  todos::add("Should I allow cores have more than one block and thus, more than one grid.");
  io::read_grid_size(&myinput, num_blocks_global, num_cells_in);
  //
  int *num_cells_per_block;
  ALLOCATE1D_INT_1ARG(num_cells_per_block, num_blocks_global);
  for (int iblock = 0; iblock < num_blocks_global; iblock++) {

    num_cells_per_block[iblock] = 1;
    for (int idir = XI; idir < DIM_MAX; idir++)
      num_cells_per_block[iblock] *= std::max( num_cells_in[iblock*DIM_MAX+idir], 1 );

  } // iblock
  int *num_cores_per_block;
  num_cores_per_block = new int[num_blocks_global];
  //
  mpi::wait_allothers("Basic grid information is initialized.");



  // assign the number of cores to each block
  decomp::assign_num_cores_crude(num_cores_per_block, num_blocks_global, num_cells_per_block);
  mpi::wait_allothers("Cores are assigned to each block.");



  // build data structure
  region = new Geometry::Generic;
  Geometry::init_region(region, &myinput, num_blocks_global, num_cells_per_block);
  io::report_geometry.initial_region(region);
  mpi::wait_allothers("region is initialized.");
  //
  // block level: block is a subset of region (i.e. region is a union of blocks)
  block = new Geometry::StructuredBlock[region->num_partitions];
  for (int iblock = 0; iblock < region->num_partitions; iblock++) {

    int iarray3[DIM_MAX];

    for (int idir = XI; idir < DIM_MAX; idir++)
      iarray3[idir] = num_cells_in[iblock*DIM_MAX+idir];

    Geometry::init_block(iblock, &block[iblock], region, &myinput, iarray3, num_cores_per_block[iblock]);

  } // iblock
  mpi::wait_allothers("block is initialized.");

  // block level: decompose blocks in each direction
  todos::add("Do a smarter and more sophisticated domain decomposition.");
  todos::add("The minimum possible number of stencil points is assumed 3: find STENCIL_MINIMUM_3.");
  int num_cores_dir[DIM_MAX];
  int num_cores_accumulated = 0;
  for (int iblock = 0; iblock < region->num_partitions; iblock++) {

    // decompose a block one-by-one
    decomp::decompDomain(myinput.how2decompDomain, myinput.file_decompDomain, iblock, num_cores_dir, block[iblock].num_cores, block[iblock].num_cells_dir);
    //
    block[iblock].num_cores = 1;
    for (int idir = XI; idir < DIM_MAX; idir++) {

      block[iblock].num_cores_dir[idir] = num_cores_dir[idir];
      block[iblock].num_cores *= num_cores_dir[idir]; // re-compute the number of cores for this block

    } // idir
    num_cores_accumulated += block[iblock].num_cores;
    block[iblock].num_partitions = block[iblock].num_cores; // update the number of partitions for this block

    // index map for decomposition within this block: this is where the core ordering is determined!
    int size_index_map = block[iblock].num_cores * DIM_MAX * 2; // beginning and ending indices (thus, 2) for each partition (thus, block[iblock].num_cores) in each direction (thus, DIM_MAX)
    int *index_map = new int[size_index_map]; 
    for (int i = 0; i < size_index_map; i++)
      index_map[i] = -1;
    Geometry::gen_index_map(&block[iblock], index_map);

    // grid level: some initialization
    if (Geometry::check_if_this_is_my_block(iblock, block)) { // if the current core belongs to this block

      grid = new Geometry::StructuredGrid;
      Geometry::init_grid(grid, iblock, block, index_map);

    } // Geometry::check_if_this_is_my_block

    DEALLOCATE_1DPTR(index_map);

  } // iblock
  //
  if ( num_cores_accumulated != mpi::nprocs )
    mpi::graceful_exit("The requested number of cores is different from what is re-computed using the provided domain decomposition.");
  io::report_parallel.initial_workLoad(grid, region->num_cells);
  io::report_geometry.initial_block_structured(block);
  io::report_geometry.initial_grid_structured(grid);
  //
  mpi::wait_allothers("block is decomposed.");



  // local communicators
  Geometry::group_cores_within_region();
  Geometry::group_cores_within_block(grid, block);
  Geometry::group_cores_head(block);
  mpi::wait_allothers("Local communicators are created.");



  // read or initialize grid points
  todos::add("Only PLOT3D file in the whole, multi-block format is supported now.");
  todos::add("Need to be able to read a PLOT3D grid file with IBLANK.");
  io::read_grid(&myinput, grid, block);



  // read or initialize state variables
  state.initialize_state(&myinput, grid);
  mpi::wait_allothers("Solution is initialized.");
  if (myinput.present_file_solution_in == TRUE) {

    io::read_solution(myinput.file_solution_in, &myinput, grid, block, &state);
    mpi::wait_allothers("The last solution is up and the simulation will restart.");

  } // myinput.present_file_solution_in
  if (myinput.present_file_mean_in == TRUE) {

    io::read_meanState(myinput.file_mean_in, &myinput, grid, block, state.num_vars_mean, state.sol_mean);
    mpi::wait_allothers("Mean state is read.");

  } // myinput.present_file_mean_in
  if (myinput.present_file_aux_in == TRUE) {

    io::read_auxvar(myinput.file_aux_in, &myinput, grid, block, state.num_vars_aux, state.sol_aux);
    mpi::wait_allothers("Auxiliary variables are read.");

  } // myinput.present_file_mean_in
  if (myinput.present_file_solution_in == TRUE ||
      myinput.present_file_mean_in == TRUE ||
      myinput.present_file_aux_in == TRUE)
    state.compute_dependent_variables(state.sol); // if the solution or the mean state is updated by reading a file, 
                                                  // compute the dependent variables once again



  // index overriding
  grid->override_local_indices(&myinput);



  // boundary
  io::read_bc(&myinput, grid, block); // read from the boundary condition file
  Geometry::additionalInit_boundary(&myinput, grid, block); // make relevant changes to the boundary data
  mpi::wait_allothers("Boundary is initialized.");



  // initialize solver
  solver::initialize(&myinput, grid, &state);
  mpi::wait_allothers("Solver is initialized.");



  // metrics
  todos::add("If purely 1-D or 2-D, derivative and filter are set to be zero in that direction; could be problematic for axisymmetric cases: find OPERATOR_IF_ONECELL.");
  todos::add("Metric computation along a plane periodic direction: there could be a huge jump in x, y, and z.");
  todos::add("Do corners cause any issues?");
  metrics::compute(&myinput, grid);
  mpi::wait_allothers("Grid metrics are computed.");


  // overset grid information
  if (myinput.do_overset == TRUE) {

    overset::initialize_overset(&myinput, grid, block);
    mpi::wait_allothers("Overset grid is initialized.");

  } // myinput.do_overset



  // initialize probe
  probe::initialize(&myinput);



  // write the initial data before time marching
  io::write_grid(&myinput, grid, block);
  io::write_metrics(&myinput, grid, block);
  io::write_solution(&myinput, grid, block, &state);
  io::write_solution_mean(&myinput, grid, block, &state);
  io::write_solution_aux(&myinput, grid, block, &state);
//  if (myinput.do_overset == TRUE) {
//
//    overset::interpolate(grid, state.sol);
//    io::write_solution(&myinput, grid, block, &state);
//    mpi::graceful_exit("Test overset interpolation is done.");
//
//  } // myinput.do_overset
  mpi::wait_allothers("Initial data are written.");



  // prepare to time march
  Geometry::clean_up_before_time_marching(grid);
  solver::precompute_something(&myinput, grid, &state);
  mpi::wait_allothers("Ready to do time-marching.");

  return;

} // initialize



void run() {

  // solve
  todos::add("Read a text file every time step to intervene the code running.");
  todos::add("Report CPU time so that the scalability can be estimated (overall/modules).");
  todos::add("Have a function pointer for computing contravariant velocities in 2-D/3-D and in time-stationary/moving coordinates.");
  todos::add("Sometimes, simulation suddenly blows up; could be some MPI issues.");
  solver::solve(&myinput, grid, block, &state);

  return;

} // run



void finalize() {

  // terminate the code
  todos::add("Additional clean-ups are needed.");
  solver::finalize();
  mpi::end_parallel_global();

  return;

} // finalize

} // simulation



namespace todos {

int num_todos = 0;

void add(std::string what_I_need_to_do) {

  std::ofstream ofs;

  num_todos++;

  if (mpi::irank == 0) {

    ofs.open(cstr_to_constchar("./report/todos"), std::ofstream::app);
    ofs << "Task #: " << std::setw(4) << num_todos << ", " << what_I_need_to_do << std::endl;
    ofs.close();

  }

  return;

} // add

} //todos
