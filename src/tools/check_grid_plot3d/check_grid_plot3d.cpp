#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "../../core/param.h"
#include "../../core/macros_inlines.h"
#include "../../io/plot3d.h"
#include "../../parallel/parallel.h"

int main(int argc, char * argv[]) {

  std::string filename;
  int num_blocks;
  int *num_cells_dir;
  std::string str_stdout;
  std::stringstream str_counter;

  // start MPI
  MPI_Init(&argc, &argv);
  mpi::init_parallel_global();

  // if nothing is specified for a command-line argument, exit
  if (argc != 2)
    mpi::graceful_exit("Specify a name for a PLOT3D grid file.");

  filename = argv[1]; // argv[0] reserved for the executable
  //std::cout << filename << std::endl;

  // print out number of blocks and number of grid points
  plot3d::stdout_grid_size(filename, PLOT3D_GRID);

  // read the file header and keept the file open
  std::ifstream ifs;
  plot3d::read_grid_header_fileKeptOpen(filename, num_blocks, num_cells_dir, ifs);

  // read the body of the file and get the min/max of each variable
  int num_vars = DIM_MAX;
  double *minmax;
  ALLOCATE1D_DOUBLE_2ARG(minmax, num_vars, 2);
  plot3d::minmax_in_gridFile(num_blocks, num_cells_dir, num_vars, minmax, ifs);

  // report
  // min/max of variables
  if (mpi::irank == 0) {
    for (int ivar = 0; ivar < num_vars; ivar++) {

      str_counter.str("");
      str_counter << std::setw(6) << ivar;
      str_stdout = "Variable " + str_counter.str() + ": minimum ";

      str_counter.str("");
      str_counter << std::setw(20) << minmax[2*ivar+MINIMUM];
      str_stdout += str_counter.str() + "; maximum ";

      str_counter.str("");
      str_counter << std::setw(20) << minmax[2*ivar+MAXIMUM];
      str_stdout += str_counter.str();

      MESSAGE_STDOUT(str_stdout);

    } // ivar
  } // mpi::irank
  MESSAGE_STDOUT("");

  // clean up
  DEALLOCATE_1DPTR(num_cells_dir);
  DEALLOCATE_1DPTR(minmax);

  mpi::end_parallel_global();

  return 0;

} // main
