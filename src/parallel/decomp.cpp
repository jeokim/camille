#include <string>
#include <cmath>

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../math/algebra.h"
#include "../io/io.h"
#include "parallel.h"
#include "decomp.h"

namespace decomp {

void assign_num_cores_crude(int *num_cores_per_block, int num_blocks_global, int *num_cells_per_block) { // DECOMP_ONEDIRECTION

  int num_cells_global;
  int cur_sum_cores;
  int iblock_maxcores, iblock_maxcells;

  if (mpi::nprocs == num_blocks_global) {

    for (int iblock = 0; iblock < num_blocks_global; iblock++)
      num_cores_per_block[iblock] = 1; // each block has to have at least one core to run
                                       // this may end up being a horrible load balance, though

    return;

  } // mpi::nprocs

  else if (mpi::nprocs < num_blocks_global)
    mpi::graceful_exit("Not ready to do this case of nprocs smaller than num_blocks_global.");

  // calculate the number of cells over the entire domain
  num_cells_global = math_algebra::sum_array_integer(num_cells_per_block, num_blocks_global);

  // assign cores according to workload estimated by number of cells
  for (int iblock = 0; iblock < num_blocks_global; iblock++) {

    num_cores_per_block[iblock] = floor( mpi::nprocs * num_cells_per_block[iblock] / num_cells_global );
    num_cores_per_block[iblock] = std::max(num_cores_per_block[iblock], 1); // there should be at least one core per block

  } // iblock

  // calculate the current sum of cores
  cur_sum_cores = math_algebra::sum_array_integer(num_cores_per_block, num_blocks_global);
  while (cur_sum_cores > mpi::nprocs) {

    iblock_maxcores = maxloc(num_cores_per_block, num_blocks_global); // a block having the largest number of cores will lose one
    num_cores_per_block[iblock_maxcores] -= 1;
    cur_sum_cores = math_algebra::sum_array_integer(num_cores_per_block, num_blocks_global);

  } // cur_sum_cores
  //
  while (cur_sum_cores < mpi::nprocs) {

    iblock_maxcells = maxloc(num_cells_per_block, num_blocks_global); // a block having the largest number of cells will gain one
    num_cores_per_block[iblock_maxcells] += 1;
    cur_sum_cores = math_algebra::sum_array_integer(num_cores_per_block, num_blocks_global);

  } // cur_sum_cores
  cur_sum_cores = math_algebra::sum_array_integer(num_cores_per_block, num_blocks_global);
  if (cur_sum_cores != mpi::nprocs)
    mpi::graceful_exit("The algorithm adjusting the number of cores failed somehow.");

  return;

} // assign_num_cores_crude



void decompDomain(std::string how2decompDomain, std::string file_decompDomain, int iblock, int num_cores_dir[DIM_MAX], int num_cores, int num_cells_dir[DIM_MAX]) {

  // this is a master function for domain decomposition
  // note that this function decompose only a single block; 
  // thus, this function should be called multiple times for a multi-block domain

  if (how2decompDomain == "DECOMP_1D")
    decomp_1dir(num_cores_dir, num_cores, num_cells_dir);

  else if (how2decompDomain == "DECOMP_FROMFILE")
    decomp_fromFile(file_decompDomain, iblock, num_cores_dir);

  else
    mpi::graceful_exit("The current domain-decomposition method is unknown.");

  return;

} // decompDomain



void decomp_1dir(int num_cores_dir[DIM_MAX], int num_cores, int num_cells_dir[DIM_MAX]) {

  int idir_max_num_cells;
  int avail_cores;
  int avg_num_cells_dir[DIM_MAX];

  for (int idir = XI; idir < DIM_MAX; idir++)
    num_cores_dir[idir] = 1;

  idir_max_num_cells = maxloc(num_cells_dir, DIM_MAX);
  avail_cores = num_cores;

  for (int idir = XI; idir < DIM_MAX; idir++)
    avg_num_cells_dir[idir] = num_cells_dir[idir];
  avg_num_cells_dir[idir_max_num_cells] = floor(num_cells_dir[idir_max_num_cells] / avail_cores);

  if (avg_num_cells_dir[idir_max_num_cells] < 3) // quit if there are less than 3 cells in the decomposing direction: STENCIL_MINIMUM_3
    mpi::graceful_exit("There should be at least 3 cells in the decomposing direction.");

  num_cores_dir[idir_max_num_cells] = std::max( avail_cores, 1 );

  if (math_algebra::product_array_integer(num_cores_dir, DIM_MAX) != num_cores)
    mpi::graceful_exit("Inconsistent number of cores after decomposing in one direction.");

  return;

} // decomp_1dir



void decomp_fromFile(std::string file_decompDomain, int iblock, int num_cores_dir[DIM_MAX]) {

  io::read_file_decompDomain(file_decompDomain, iblock, num_cores_dir);

  return;

} // decomp_fromFile

} // decomp
