#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "geometry.h"

namespace Geometry {

double *xyz_block;



void init_region(Generic *region, UserInput *input, int num_blocks_global, int *num_cells_per_block) {

  region->num_dim = input->num_dim;

  region->num_ourkinds = 1; // there is only one region
  region->num_partitions = num_blocks_global;
  region->num_cores = mpi::nprocs;

  region->id_global = 0; // there is only one region
  region->id_local = 0; // there is only one region

  region->id_parent = NONE; // no parent since region is highest in the hierarchy

  region->num_cells = math_algebra::sum_array_integer(num_cells_per_block, num_blocks_global);
  //region->num_ocells is not initialized since region does not know anything about ghost cells (only grid knows)

  return;

} // init_region



void init_block(int iblock, StructuredBlock *block, Generic *region, UserInput *input, int num_cells_dir_in_this_block[DIM_MAX], int num_cores_in_this_block) {

  block->num_dim = region->num_dim;

  block->num_ourkinds = region->num_partitions; // how many blocks there are within the upper-level data structure
  block->num_partitions = num_cores_in_this_block;
  block->num_cores = num_cores_in_this_block;

  block->id_global = iblock;
  block->id_local = iblock;

  block->id_parent = 0; // all blocks has the same parent region

  block->num_cells = 1;
  for (int idir = XI; idir < DIM_MAX; idir++)
    block->num_cells *= std::max( num_cells_dir_in_this_block[idir], 1 );
  //block->num_ocells is not initialized since block does not know anything about ghost cells (only grid knows)

  //

  block->num_cells_ghost = input->num_cells_ghost;

  for (int idir = XI; idir < DIM_MAX; idir++)
    block->periodic[idir] = NONPERIODIC;

  for (int idir = XI; idir < DIM_MAX; idir++) {

    block->num_cells_dir[idir] = num_cells_dir_in_this_block[idir];
    //block->num_ocells_dir[] is not initialized since block does not know anything about ghost cells (only grid knows)
    //
    block->is[idir] = 0;
    block->ie[idir] = num_cells_dir_in_this_block[idir] - 1;
    //block->iso[] and block->ieo[] are not initialized since block does not know anything about ghost cells (only grid knows)

  } // idir

  block->init_idx1D(block->num_cells_dir[XI], block->num_cells_dir[ETA], block->num_cells_dir[ZETA]);

  return;

} // init_block



void init_grid(StructuredGrid *mygrid, int cur_block, StructuredBlock *block, int *index_map) {

  int num_cores_sofar = 0;
  for (int iblock = 0; iblock < cur_block; iblock++)
    num_cores_sofar += block[iblock].num_cores;

  mygrid->num_dim = block[cur_block].num_dim;

  mygrid->num_ourkinds = block[cur_block].num_partitions;
  mygrid->num_partitions = 1; // grid is a minimum element and thus indivisible
  mygrid->num_cores = 1; // grid is a minimum element and thus indivisible

  mygrid->id_global = mpi::irank;
  mygrid->id_local = mpi::irank - num_cores_sofar;

  mygrid->id_parent = cur_block;

  // though block does not use any ghost cells, it is good to keep ghost cell information at the block level and distribute to its grids so that all the grids within a block share the same information
  mygrid->num_cells_ghost = block[cur_block].num_cells_ghost;

  for (int idir = XI; idir < DIM_MAX; idir++)
    mygrid->periodic[idir] = block->periodic[idir];

  // generate index maps for grid (parent/local)
  for (int idir = XI; idir < DIM_MAX; idir++) {

    // task 1: get indices at the "parent-block" level of this "grid"
    // get index_map as is, which contains non-overlapping indices for "block" (not "grid") without ghost cells
    mygrid->is_in_parent[idir] = index_map[mygrid->id_local*DIM_MAX*2 + idir*2 + 0];
    mygrid->ie_in_parent[idir] = index_map[mygrid->id_local*DIM_MAX*2 + idir*2 + 1];

    if (mygrid->is_in_parent[idir] == mygrid->ie_in_parent[idir]) { // for purely 1-D or 2-D cases

      if (mygrid->is_in_parent[idir] != 0)
        mpi::graceful_exit("Grid index for a 2-D configuration is inconsistent");

      // include ghost cells (zero in this case)
      mygrid->iso_in_parent[idir] = mygrid->is_in_parent[idir];
      mygrid->ieo_in_parent[idir] = mygrid->ie_in_parent[idir];

    }
    else {

      // pad ghost cells only to internal boundaries doing message passing; 
      // if a boundary is a block's boundary (subject to boundary condition or overset-grid interpolation), no ghost cells are added (i.e. iso = is or ieo = ie)
      mygrid->iso_in_parent[idir] = std::max( mygrid->is_in_parent[idir] - mygrid->num_cells_ghost, block[cur_block].is[idir]);
      mygrid->ieo_in_parent[idir] = std::min( mygrid->ie_in_parent[idir] + mygrid->num_cells_ghost, block[cur_block].ie[idir]);

    } // mygrid->is_in_parent[idir]

    // task 2: get the number of cells: same for both parent-block and local levels
    mygrid->num_cells_dir[idir] = mygrid->ie_in_parent[idir] - mygrid->is_in_parent[idir] + 1;
    mygrid->num_ocells_dir[idir] = mygrid->ieo_in_parent[idir] - mygrid->iso_in_parent[idir] + 1;

    // task 3: get indices at the local level of this "grid"
    // start from where ghost cells are already included
    mygrid->iso[idir] = 0;
    mygrid->ieo[idir] = mygrid->num_ocells_dir[idir] - 1;

    if (mygrid->iso[idir] == mygrid->ieo[idir]) { // for purely 1-D or 2-D cases

      if (mygrid->iso[idir] != 0)
        mpi::graceful_exit("Grid index for a 2-D configuration is inconsistent");

      // get the physical point indices
      mygrid->is[idir] = mygrid->iso[idir];
      mygrid->ie[idir] = mygrid->ieo[idir];

    }
    else {

      // get the physical point indices
      // pad ghost cells only to internal boundaries doing message passing; 
      // if a boundary is a block's boundary (subject to boundary condition or overset-grid interpolation), no ghost cells are added (i.e. iso = is or ieo = ie)
      mygrid->is[idir] = mygrid->iso[idir];
      if (mygrid->is_in_parent[idir] != block[cur_block].is[idir])
        mygrid->is[idir] += mygrid->num_cells_ghost;
      //
      mygrid->ie[idir] = mygrid->ieo[idir];
      if (mygrid->ie_in_parent[idir] != block[cur_block].ie[idir])
        mygrid->ie[idir] -= mygrid->num_cells_ghost;

    } // mygrid->is_in_parent[idir]

    //// task 4: repeat above steps for nodes
    //if (mygrid->num_cells_dir[idir] == 1) { // for purely 1-D or 2-D cases
    //
    //  mygrid->num_nodes_dir[idir] = 1;
    //  mygrid->num_onodes_dir[idir] = 1; // since it is 2-D, no ghost cell in the third direction
    //
    //}
    //else {
    //
    //  mygrid->num_nodes_dir[idir] = mygrid->num_cells_dir[idir] + 1;
    //  mygrid->num_onodes_dir[idir] = mygrid->num_ocells_dir[idir] + 1;
    //
    //} // mygrid->num_cells_dir[idir]
    //
    //// local indices
    //mygrid->js[idir] = mygrid->is[idir];
    //mygrid->je[idir] = mygrid->js[idir] + mygrid->num_nodes_dir[idir] - 1;
    //
    //mygrid->jso[idir] = mygrid->iso[idir];
    //mygrid->jeo[idir] = mygrid->jso[idir] + mygrid->num_onodes_dir[idir] - 1;
    //
    //// indices in its parent block
    //mygrid->js_in_parent[idir] = mygrid->is_in_parent[idir];
    //mygrid->je_in_parent[idir] = mygrid->js_in_parent[idir] + mygrid->num_nodes_dir[idir] - 1;
    //
    //mygrid->jso_in_parent[idir] = mygrid->iso_in_parent[idir];
    //mygrid->jeo_in_parent[idir] = mygrid->jso_in_parent[idir] + mygrid->num_onodes_dir[idir] - 1;

  } // idir

  mygrid->num_cells = mygrid->num_cells_dir[XI] * mygrid->num_cells_dir[ETA] * mygrid->num_cells_dir[ZETA];
  //mygrid->num_nodes = mygrid->num_nodes_dir[XI] * mygrid->num_nodes_dir[ETA] * mygrid->num_nodes_dir[ZETA];

  mygrid->num_ocells = mygrid->num_ocells_dir[XI] * mygrid->num_ocells_dir[ETA] * mygrid->num_ocells_dir[ZETA];
  //mygrid->num_onodes = mygrid->num_onodes_dir[XI] * mygrid->num_onodes_dir[ETA] * mygrid->num_onodes_dir[ZETA];

  for (int idir = XI; idir < DIM_MAX; idir++)
    mygrid->num_cores_dir[idir] = 1; // grid is a minimum element and thus indivisible

  mygrid->init_idx1D(mygrid->num_ocells_dir[XI], mygrid->num_ocells_dir[ETA], mygrid->num_ocells_dir[ZETA]);

  mygrid->cell = new GridPoint[mygrid->num_ocells];

  init_iblank(mygrid);

  return;

} // init_grid



void init_iblank(StructuredGrid *mygrid) {

  for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
    mygrid->cell[l0].iblank = NOTBLANKED; // by default, all grid points are fluid cells

  return;

} // init_grid



void gen_index_map(StructuredBlock *block, int *index_map) {

  int avg_num_cells, rem_num_cells;
  int is, ie;
  int counter_core;

  for (int idir = XI; idir < DIM_MAX; idir++) {

    avg_num_cells = floor( block->num_cells_dir[idir] / block->num_cores_dir[idir] );
    rem_num_cells = block->num_cells_dir[idir] - avg_num_cells * block->num_cores_dir[idir]; // remaining number of cells after evenly distributing cells to all cores in this direction; 
                                                                                             // should be smaller than the number of cores in this direction
    ie = -1;
    for (int icore = 0; icore < block->num_cores_dir[idir]; icore++) {

      is = std::max(ie + 1, 0);
      ie = std::max(is + avg_num_cells - 1, 0);

      if (rem_num_cells > 0) { // each core gets one additional cell if there are any remaining cells

        ie += 1;
        rem_num_cells -= 1;

      } // rem_num_cells

      //std::cout << "Block: " << block->id_global << ", idir: " << idir << ", icore: " << icore << ", is, ie: " << is << ", " << ie << std::endl;

      counter_core = 0;
      switch ( idir ) {
      case XI:
        for (int k = 0; k < block->num_cores_dir[ZETA]; k++)
          for (int j = 0; j < block->num_cores_dir[ETA]; j++)
            for (int i = 0; i < block->num_cores_dir[XI]; i++) {

              if (icore == i) {

                index_map[counter_core*DIM_MAX*2 + XI*2 + 0] = is;
                index_map[counter_core*DIM_MAX*2 + XI*2 + 1] = ie;

              } // icore

              counter_core++;

            } // i

        break;

      case ETA:
        for (int k = 0; k < block->num_cores_dir[ZETA]; k++)
          for (int j = 0; j < block->num_cores_dir[ETA]; j++)
            for (int i = 0; i < block->num_cores_dir[XI]; i++) {

              if (icore == j) {

                index_map[counter_core*DIM_MAX*2 + ETA*2 + 0] = is;
                index_map[counter_core*DIM_MAX*2 + ETA*2 + 1] = ie;

              } // icore

              counter_core++;

            } // i

        break;

      case ZETA:
        for (int k = 0; k < block->num_cores_dir[ZETA]; k++)
          for (int j = 0; j < block->num_cores_dir[ETA]; j++)
            for (int i = 0; i < block->num_cores_dir[XI]; i++) {

              if (icore == k) {

                index_map[counter_core*DIM_MAX*2 + ZETA*2 + 0] = is;
                index_map[counter_core*DIM_MAX*2 + ZETA*2 + 1] = ie;

              } // icore

              counter_core++;

            } // i

        break;

      } // idir
    } // icore
  } // idir

  if (counter_core != block->num_cores)
    mpi::graceful_exit("Number of cores for block decomposition is inconsistent.");

  return;

} // gen_index_map



int check_if_this_is_my_block(int cur_block, StructuredBlock *block) {

  int num_cores_sofar;
  int it_is;

  num_cores_sofar = 0;
  it_is = FALSE;

  for (int iblock = 0; iblock < cur_block; iblock++)
    num_cores_sofar += block[iblock].num_cores;

  if (mpi::irank >= num_cores_sofar && mpi::irank < num_cores_sofar + block[cur_block].num_cores)
    it_is = TRUE;

  return it_is;

} // check_if_this_is_my_block



void clean_up_before_time_marching(StructuredGrid *mygrid) {

  DEALLOCATE_1DPTR(xyz_block);

  // deallocate the temporary in/out boundaries of a buffer zone
  int num_bufferZones = mygrid->num_bufferZones;
  if (num_bufferZones > 0) {
    for (int ibuffer = FIRST; ibuffer < num_bufferZones; ibuffer++) {

      StructuredBufferZone *bufferZone_cur = &(mygrid->bufferZone[ibuffer]);

      for (int iside = IN; iside <= OUT; iside++) {

        StructuredBoundary *boundary_cur = &(bufferZone_cur->boundary[iside]);

        for (int idir = 0; idir < DIM_MAX - 1; idir++)
          DEALLOCATE_1DPTR(boundary_cur->xyz[idir]);

      } // iside
    } // ibuffer
  } // num_bufferZones

  return;

} // clean_up_before_time_marching



void fill_in_xyz_block(double *xyz_block_in, int data_size) {

  xyz_block = new double[data_size];
  for (int i = 0; i < data_size; i++)
    xyz_block[i] = xyz_block_in[i];

  return;

} // fill_in_xyz_block



void StructuredBlock::init_idx1D(int nXi_in, int nEta_in, int nZeta_in) {

  nXi = nXi_in;
  nEta = nEta_in;
  nZeta = nZeta_in;

  return;

} // StructuredBlock::init_idx1D



void StructuredGrid::hardwire_gridPoint(UserInput *input, StructuredBlock *myblock) {

  // manually specify coordinate locations at physical points for each block

//  if (this->id_parent == 0) {
//
//    double r_in = 0.1;
//    double r_out = 1.0;
//    double dr = (r_out - r_in) / static_cast<double>(myblock->num_cells_dir[XI]); // xi corresponds to a radial direction
//    r_in += dr * 0.5;
//    double theta_0 = 0.0;
//    double theta_1 = math_constants::pi / 4.0;
//    double dtheta = fabs(theta_1 - theta_0) / static_cast<double>(myblock->num_cells_dir[ETA]); // eta corresponds to an azimuthal direction
//
//    for (int k = this->iso[ZETA]; k < this->ieo[ZETA] + 1; k++) {
//      int k_in_parent = this->iso_in_parent[ZETA] + k;
//
//      for (int j = this->iso[ETA]; j < this->ieo[ETA] + 1; j++) {
//        int j_in_parent = this->iso_in_parent[ETA] + j;
//
//        for (int i = this->iso[XI]; i < this->ieo[XI] + 1; i++) {
//          int i_in_parent = this->iso_in_parent[XI] + i;
//
//          int l0 = this->idx1D(i, j, k);
//
//          double r = r_in + static_cast<double>(i_in_parent) * dr;
//          double theta = theta_0 + static_cast<double>(j_in_parent) * dtheta;
//
//          r *= 1.0 + (theta - theta_0) * 0.5; // radius increases as the eta-index increases
//          this->cell[l0].xyz[XDIR] = r * cos(theta);
//          this->cell[l0].xyz[YDIR] = r * sin(theta);
//          this->cell[l0].xyz[ZDIR] = 0.0;
//
//        } // i
//      } // j
//    } // k
//  }

  double xL = 2.0, yL = 0.5;
  double dx = xL / static_cast<double>(myblock->num_cells_dir[XI]), dy = yL / static_cast<double>(myblock->num_cells_dir[ETA]);
  double x0 = -xL / 2.0 + dx / 2.0, y0 = -yL / 2.0 + dy / 2.0;

  // manually specify coordinate locations at physical points for each block
  if (this->id_parent == 0) {

    for (int k = this->iso[ZETA]; k < this->ieo[ZETA] + 1; k++) {
      int k_in_parent = this->iso_in_parent[ZETA] + k;

      for (int j = this->iso[ETA]; j < this->ieo[ETA] + 1; j++) {
        int j_in_parent = this->iso_in_parent[ETA] + j;

        for (int i = this->iso[XI]; i < this->ieo[XI] + 1; i++) {
          int i_in_parent = this->iso_in_parent[XI] + i;

          int l0 = this->idx1D(i, j, k);

          this->cell[l0].xyz[XDIR] = x0 + static_cast<double>(i_in_parent) * dx;
          this->cell[l0].xyz[YDIR] = y0 + static_cast<double>(j_in_parent) * dy;

        } // i
      } // j
    } // k
  } // this->id_parent
  else
  if (this->id_parent == 1) {

    for (int k = this->iso[ZETA]; k < this->ieo[ZETA] + 1; k++) {
      int k_in_parent = this->iso_in_parent[ZETA] + k;

      for (int j = this->iso[ETA]; j < this->ieo[ETA] + 1; j++) {
        int j_in_parent = this->iso_in_parent[ETA] + j;

        for (int i = this->iso[XI]; i < this->ieo[XI] + 1; i++) {
          int i_in_parent = this->iso_in_parent[XI] + i;

          int l0 = this->idx1D(i, j, k);

          this->cell[l0].xyz[XDIR] = (x0 + 0.0) + static_cast<double>(i_in_parent) * dx;
          this->cell[l0].xyz[YDIR] = (y0 + 0.7) + static_cast<double>(j_in_parent) * dy;
          this->cell[l0].xyz[ZDIR] = 0.0;

        } // i
      } // j
    } // k

  } // this->id_parent

  return;

} // StructuredGrid::hardwire_gridPoint



void StructuredBoundary::init_idx1D(int nXi_in, int nEta_in, int nZeta_in) {

  nXi = nXi_in;
  nEta = nEta_in;
  nZeta = nZeta_in;

  return;

} // StructuredBoundary::init_idx1D



void StructuredGrid::initialize_boundaryCondition(int num_boundaryCondition_nonperiodic_in) {

  this->num_boundaryCondition_nonperiodic = num_boundaryCondition_nonperiodic_in;

  if (num_boundaryCondition_nonperiodic_in > 0)
    this->boundaryCondition = new StructuredBoundaryCondition[this->num_boundaryCondition_nonperiodic];

  return;

} // StructuredGrid::initialize_boundaryCondition



void StructuredGrid::initialize_bufferZone(int num_bufferZones_in) {

  this->num_bufferZones = num_bufferZones_in;

  if (num_bufferZones_in > 0)
    this->bufferZone = new StructuredBufferZone[this->num_bufferZones];

  return;

} // StructuredGrid::initialize_bufferZone



int StructuredGrid::check_if_this_is_my_point(int num_dim, double xyz_in[DIM_MAX], int *&ijk) {

  // given a point, determine if this grid contains it and return the enclosing cell's ijk indices at the grid level

  double vec0[DIM_MAX], vec1[DIM_MAX];
  double eps_smallAngle = 0.05; //0.0; // threshold value for how small sin(angle[rad]) is small
if (mpi::irank == 1 || mpi::irank == 2)
std::cout << "I came in." << std::endl;
  // Test 1: if this point is out of bound
  //         a quick way to reject points
  int out_of_bound = FALSE;
  double coordinate_min, coordinate_max;
  double *tmp;
  ALLOCATE1D_DOUBLE_1ARG(tmp,this->num_cells);
  for (int idir = XDIR; idir < DIM_MAX; idir++) {
    int count = 0;
    for (int k = this->iso[ZETA]; k <= this->ieo[ZETA]; k++)
      for (int j = this->iso[ETA]; j <= this->ieo[ETA]; j++)
        for (int i = this->iso[XI]; i <= this->ieo[XI]; i++) {

          int l0 = this->idx1D(i, j, k);

          tmp[count++] = this->cell[l0].xyz[idir]; // dump either x, y, or z into a 1-D array

        } // i
    assert(count == this->num_ocells);

    coordinate_min = math_algebra::minval(tmp,this->num_ocells); // minimum of either x, y, or z
    coordinate_max = math_algebra::maxval(tmp,this->num_ocells); // maximum of either x, y, or z

    if ((xyz_in[idir]-coordinate_min)*(xyz_in[idir]-coordinate_max) > 0) {
      out_of_bound = TRUE;
      break;
    } // (xyz_in[idir]-coordinate_min)*(xyz_in[idir]-coordinate_max)
  } // idir
  DEALLOCATE_1DPTR(tmp);
  if (out_of_bound == TRUE)
    return FALSE;
if (mpi::irank == 1 || mpi::irank == 2)
std::cout << "Test 1." << std::endl;
  // Test 2: at least, this point is in-bound; need a more detailed search
  //         go over every cell made out of 4 (in 2D) or 8 (in 3D) points
  int found = FALSE;
  if (num_dim == 2) {
    int num_vertices_per_cv = 4; // 2-D rectangle cell (or control volume) has 4 vertices
    int vertex_list[num_vertices_per_cv];

    int k = this->is[ZETA];
    for (int j = this->iso[ETA]; j <= this->ieo[ETA]-1; j++) // -1 since searching over cells, not points
      for (int i = this->iso[XI]; i <= this->ieo[XI]-1; i++) { // -1 since searching over cells, not points

        if (found == TRUE)
          break;

        // vertices are ordered in the counter clockwise direction
        // so if you walk along a cell's periphery, an interior point should be always on your left
        vertex_list[FIRST]  = this->idx1D(i,   j,   k);
        vertex_list[SECOND] = this->idx1D(i+1, j,   k);
        vertex_list[THIRD]  = this->idx1D(i+1, j+1, k);
        vertex_list[FOURTH] = this->idx1D(i,   j+1, k);

        // if any of points is not a fluid point (that is, blanked or interpolation point), reject the corresponding cell
        int cell_is_fluid = TRUE;
        for (int ivertex = FIRST; ivertex < num_vertices_per_cv; ivertex++) {
          int l0 = vertex_list[ivertex];
          if (this->cell[l0].iblank != NOTBLANKED)
            cell_is_fluid = FALSE;
        } // ivertex
        if (cell_is_fluid == FALSE)
          continue;

        // actual search
        int point_inside = TRUE;
        int ivertex = FIRST;
        while ( point_inside && ivertex < num_vertices_per_cv ) {

          // get the indices of two neighboring vertices
          int ivertex0 = ivertex;
          int ivertex1 = (ivertex0 + 1) % num_vertices_per_cv;

          // get their 1-D indices in grid
          int l0 = vertex_list[ivertex0];
          int l1 = vertex_list[ivertex1];

          // form a pair of vectors, vec0 and vec1
          for (int idim = XDIR; idim < num_dim; idim++) {

            vec0[idim] = this->cell[l1].xyz[idim] - this->cell[l0].xyz[idim];
            vec1[idim] = xyz_in[idim]             - this->cell[l0].xyz[idim];

          } // idim

          // note that this cross_product is equivalent to a sine of the angles between the two vectors
          double cross_product = math_matrix::cross_product_normalized(vec0, vec1, num_dim);

          if (cross_product < -eps_smallAngle) // the inquiry point lies outside of this cell
            point_inside = FALSE;

          else if (cross_product >= -eps_smallAngle) // the inquiry point is still inside; keep searching
            point_inside = TRUE;

          ivertex++;

        } // point_inside

        if ( point_inside ) { // a cell enclosing our inquiry point is found
          found = TRUE;
          ijk[XI] = i; ijk[ETA] = j; ijk[ZETA] = k;
        } // point_inside
      } // i
if (mpi::irank == 1 || mpi::irank == 2)
std::cout << "Test 2." << std::endl;

  } // num_dim
  else if (num_dim == 3) {

    mpi::graceful_exit("check_if_this_is_my_point is not implemented for 3D yet.");

  } // num_dim

  return found;

} // StructuredGrid::check_if_this_is_my_point



void StructuredBoundary::initialize(int num_dim_in, int is_in_parent_in[], int ie_in_parent_in[], int idir_boundary, int iend_boundary) {

  this->num_dim = num_dim_in;

  this->num_cells_dir[idir_boundary] = 1;
  for (int idir = 0; idir < DIM_MAX - 1; idir++) {

    int idir_other = dir_other[idir_boundary][idir];
    this->num_cells_dir[idir_other] = ie_in_parent_in[idir_other] - is_in_parent_in[idir_other] + 1;

  } // idir
//  std::cout << "Rank: " << mpi::irank << ", idir_boundary: " << idir_boundary << 
//            ", # of cells: " << this->num_cells_dir[XI] << " " 
//                             << this->num_cells_dir[ETA] << " " 
//                             << this->num_cells_dir[ZETA] << " " << std::endl;

  this->num_cells = math_algebra::product_array_integer(this->num_cells_dir, DIM_MAX);

  // my local indices
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is[idir] = 0;
    this->ie[idir] = this->num_cells_dir[idir] - 1;

  } // idir
//  std::cout << "Rank: " << mpi::irank << ", idir_boundary: " << idir_boundary << ", iend_boundary: " << iend_boundary
//            << ", indices: " << this->is[XI] << " " << this->ie[XI] << ", " 
//                             << this->is[ETA] << " " << this->ie[ETA] << ", " 
//                             << this->is[ZETA] << " " << this->ie[ZETA] << std::endl;

  // indices at its parent grid level (for a direct access)
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is_in_parent[idir] = is_in_parent_in[idir];
    this->ie_in_parent[idir] = ie_in_parent_in[idir];

  } // idir
  if (iend_boundary == LEFT) {

    this->is_in_parent[idir_boundary] = is_in_parent_in[idir_boundary];
    this->ie_in_parent[idir_boundary] = is_in_parent_in[idir_boundary];

  } // iend_boundary
  if (iend_boundary == RIGHT) {

    this->is_in_parent[idir_boundary] = ie_in_parent_in[idir_boundary];
    this->ie_in_parent[idir_boundary] = ie_in_parent_in[idir_boundary];

  } // iend_boundary
//  std::cout << "Rank: " << mpi::irank << ", idir_boundary: " << idir_boundary << ", iend_boundary: " << iend_boundary
//            << ", indices: " << this->is_in_parent[XI] << " " << this->ie_in_parent[XI] << ", " 
//                             << this->is_in_parent[ETA] << " " << this->ie_in_parent[ETA] << ", " 
//                             << this->is_in_parent[ZETA] << " " << this->ie_in_parent[ZETA] << std::endl;

  this->init_idx1D(this->num_cells_dir[XI], this->num_cells_dir[ETA], this->num_cells_dir[ZETA]);

  for (int idir = XI; idir < DIM_MAX; idir++)
    this->xyz[idir] = new double[this->num_cells];

  return;

} // StructuredBoundary::initialize



void StructuredBoundaryCondition::initialize(int num_dim_in, int iso_in_parent[], int ieo_in_parent[], int idir_boundary, int iend_boundary, int itype_boundary) {

  this->num_dim = num_dim_in;

  this->num_cells_dir[idir_boundary] = 1;
  for (int idir = 0; idir < DIM_MAX - 1; idir++) {

    int idir_other = dir_other[idir_boundary][idir];
    this->num_cells_dir[idir_other] = ieo_in_parent[idir_other] - iso_in_parent[idir_other] + 1;

  } // idir
//  std::cout << "Rank: " << mpi::irank << ", local rank: " << mpi::irank_block << ", idir_boundary: " << idir_boundary << 
//            ", # of cells: " << this->num_cells_dir[XI] << " " 
//                             << this->num_cells_dir[ETA] << " " 
//                             << this->num_cells_dir[ZETA] << " " << std::endl;

  this->num_cells = math_algebra::product_array_integer(this->num_cells_dir, DIM_MAX);

  // my local indices
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is[idir] = 0;
    this->ie[idir] = this->num_cells_dir[idir] - 1;

  } // idir
//  std::cout << "Rank: " << mpi::irank << ", local rank: " << mpi::irank_block << ", idir_boundary: " << idir_boundary << ", iend_boundary: " << iend_boundary
//            << ", indices: " << this->is[XI] << " " << this->ie[XI] << ", " 
//                             << this->is[ETA] << " " << this->ie[ETA] << ", " 
//                             << this->is[ZETA] << " " << this->ie[ZETA] << std::endl;

  // indices at its parent grid level (for a direct access)
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is_in_parent[idir] = iso_in_parent[idir];
    this->ie_in_parent[idir] = ieo_in_parent[idir];

  } // idir
  if (iend_boundary == LEFT) {

    this->is_in_parent[idir_boundary] = iso_in_parent[idir_boundary];
    this->ie_in_parent[idir_boundary] = iso_in_parent[idir_boundary];

  } // iend_boundary
  if (iend_boundary == RIGHT) {

    this->is_in_parent[idir_boundary] = ieo_in_parent[idir_boundary];
    this->ie_in_parent[idir_boundary] = ieo_in_parent[idir_boundary];

  } // iend_boundary
//  std::cout << "Rank: " << mpi::irank << ", local rank: " << mpi::irank_block << ", idir_boundary: " << idir_boundary << ", iend_boundary: " << iend_boundary
//            << ", indices: " << this->is_in_parent[XI] << " " << this->ie_in_parent[XI] << ", " 
//                             << this->is_in_parent[ETA] << " " << this->ie_in_parent[ETA] << ", " 
//                             << this->is_in_parent[ZETA] << " " << this->ie_in_parent[ZETA] << std::endl;

  this->init_idx1D(this->num_cells_dir[XI], this->num_cells_dir[ETA], this->num_cells_dir[ZETA]);

  set_type(idir_boundary, iend_boundary, BOUNDARY_BC, itype_boundary);

  return;

} // StructuredBoundaryCondition::initialize



void StructuredBoundaryCondition::set_type(int which_dir_in, int which_end_in, int which_boundary_in, int which_model_in) {

  this->which_dir = which_dir_in;
  this->which_end = which_end_in;
  this->which_boundary = which_boundary_in;
  this->which_model = which_model_in;

  return;

} // StructuredBoundaryCondition::set_type



void StructuredBufferZone::initialize(int num_dim_in, int iso_in_parent[], int ieo_in_parent[], int idir_boundary, int iend_boundary, int itype_buffer) {

  this->num_dim = num_dim_in;

  for (int idir = XI; idir < DIM_MAX; idir++)
    this->num_cells_dir[idir] = ieo_in_parent[idir] - iso_in_parent[idir] + 1;

//  std::cout << "Rank: " << mpi::irank << ", local rank: " << mpi::irank_block << ", idir_boundary: " << idir_boundary << 
//            ", # of cells: " << this->num_cells_dir[XI] << " " 
//                             << this->num_cells_dir[ETA] << " " 
//                             << this->num_cells_dir[ZETA] << " " << std::endl;

  this->num_cells = math_algebra::product_array_integer(this->num_cells_dir, DIM_MAX);

  // my local indices
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is[idir] = 0;
    this->ie[idir] = this->num_cells_dir[idir] - 1;

  } // idir
//  std::cout << "Rank: " << mpi::irank << ", local rank: " << mpi::irank_block << ", idir_boundary: " << idir_boundary << ", iend_boundary: " << iend_boundary
//            << ", indices: " << this->is[XI] << " " << this->ie[XI] << ", " 
//                             << this->is[ETA] << " " << this->ie[ETA] << ", " 
//                             << this->is[ZETA] << " " << this->ie[ZETA] << std::endl;

  // indices at its parent grid level (for a direct access)
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is_in_parent[idir] = iso_in_parent[idir];
    this->ie_in_parent[idir] = ieo_in_parent[idir];

  } // idir
//  std::cout << "Rank: " << mpi::irank << ", local rank: " << mpi::irank_block << ", idir_boundary: " << idir_boundary << ", iend_boundary: " << iend_boundary
//            << ", indices: " << this->is_in_parent[XI] << " " << this->ie_in_parent[XI] << ", " 
//                             << this->is_in_parent[ETA] << " " << this->ie_in_parent[ETA] << ", " 
//                             << this->is_in_parent[ZETA] << " " << this->ie_in_parent[ZETA] << std::endl;

  this->init_idx1D(this->num_cells_dir[XI], this->num_cells_dir[ETA], this->num_cells_dir[ZETA]);

  set_type(idir_boundary, iend_boundary, BOUNDARY_BUFFERZONE, itype_buffer);

  return;

} // StructuredBufferZone::initialize



void StructuredBufferZone::set_parameters_4bufferZone(UserInput *myinput, int is_of_bufferZone[], int ie_of_bufferZone[]) {

  this->buffer_polynomial_order = myinput->buffer_polynomial_order;
  this->buffer_constant = myinput->buffer_constant;

  this->boundary = new StructuredBoundary[2]; // left and right boundaries of the buffer zone
                                              // not just my portion of the buffer zone, but before the zone is decomposed
  switch ( this->which_end ) {
  case LEFT: // if a buffer zone is on the left-hand-side, disturbances flow IN from its RIGHT and flow OUT on its LEFT
    this->boundary[IN ].initialize(this->num_dim, is_of_bufferZone, ie_of_bufferZone, this->which_dir, RIGHT);
    this->boundary[OUT].initialize(this->num_dim, is_of_bufferZone, ie_of_bufferZone, this->which_dir, LEFT );

    break;

  case RIGHT: // if a buffer zone is on the right-hand-side, disturbances flow IN from its LEFT and flow OUT on its RIGHT
    this->boundary[IN ].initialize(this->num_dim, is_of_bufferZone, ie_of_bufferZone, this->which_dir, LEFT );
    this->boundary[OUT].initialize(this->num_dim, is_of_bufferZone, ie_of_bufferZone, this->which_dir, RIGHT);

    break;

  } // this->which_end

  this->buffer_strength = new double[this->num_cells];
  for (int lb = 0; lb < this->num_cells; lb++)
    (this->buffer_strength)[lb] = 0.0;

  return;

} // StructuredBufferZone::set_parameters_4bufferZone



void group_cores_within_region() {

  mpi::comm_region = MPI_COMM_WORLD; // equivalent to MPI_COMM_WORLD

  MPI_Comm_group(mpi::comm_region, &mpi::group_region);
  MPI_Group_rank(mpi::group_region, &mpi::irank_region); // equivalent to irank
  MPI_Group_size(mpi::group_region, &mpi::nprocs_region); // equivalent to nprocs

  mpi::wait_allothers();

  //std::cout << "Rank: " << mpi::irank << ", irank_region: " << mpi::irank_region << ", nprocs_region: " << mpi::nprocs_region << std::endl;

  return;

} // group_cores_within_region



void group_cores_within_block(StructuredGrid *mygrid, StructuredBlock *block) {

  int head_rank = mpi::irank - mygrid->id_local; // head rank of my block
  int n = mygrid->num_ourkinds;
  int *ranks = new int[n];

  for (int i = 0; i < n; i++)
    ranks[i] = head_rank + i;

  MPI_Group_incl(mpi::group_region, n, ranks, &mpi::group_block);
  MPI_Group_rank(mpi::group_block, &mpi::irank_block);
  MPI_Group_size(mpi::group_block, &mpi::nprocs_block);

  MPI_Comm_create(mpi::comm_region, mpi::group_block, &mpi::comm_block);

  DEALLOCATE_1DPTR(ranks);

  mpi::wait_allothers();

  //std::cout << "Rank: " << mpi::irank << ", irank_block: " << mpi::irank_block << ", nprocs_block: " << mpi::nprocs_block << std::endl;

  //{
  //  int send_message = 4 * mpi::irank + 1136; // some arbitrary number in a function of a global rank
  //  int tmp = -mpi::irank;
  //  for (int id_local = 1; id_local < mygrid->num_ourkinds; id_local++) { // starting from 1 to exclude a head node
  //    if (mygrid->id_local == id_local) { // if not a head node, send to a head node and wait
  //        //std::cout << "Rank: " << mpi::irank << ", id_local: " << id_local << std::endl;
  //        MPI_Send(&send_message, 1, MPI_INT, 0, id_local, mpi::comm_block); // 0 means a head node in this group
  //    }
  //    else if (mygrid->id_local == 0) { // if I am a head node, always receive and print
  //      //std::cout << "Rank: " << mpi::irank << ", message (before): " << tmp << " from rank "<< head_rank + id_local << std::endl;
  //      MPI_Recv(&tmp, 1, MPI_INT, id_local, id_local, mpi::comm_block, mpi::status);
  //      //std::cout << "Rank: " << mpi::irank << ", message (after): " << tmp << " from rank "<< head_rank + id_local << std::endl;
  //    } // mygrid->id_local
  //  } // id_local
  //  std::cout << "Rank: " << mpi::irank << ", tmp: " << tmp << std::endl;
  //}

  label_cores_within_block(mygrid, &block[mygrid->id_parent]);

  return;

} // group_cores_within_block



void group_cores_head(StructuredBlock *block) {

  // only group cores with head ranks

  int n = block[0].num_ourkinds; // the number of head ranks is equal to the number of blocks
  int *ranks = new int[n];

  int num_blocks = n;
  int num_cores_sofar = 0;
  for (int iblock = 0; iblock < num_blocks; iblock++) {

    ranks[iblock] = num_cores_sofar;
    num_cores_sofar += block[iblock].num_cores;

  } // iblock

  MPI_Group_incl(mpi::group_region, n, ranks, &mpi::group_head);
  MPI_Group_rank(mpi::group_head, &mpi::irank_head);
  MPI_Group_size(mpi::group_head, &mpi::nprocs_head);

  MPI_Comm_create(mpi::comm_region, mpi::group_head, &mpi::comm_head);

//  if (mpi::irank_block == 0)
//    std::cout << "Rank: " << mpi::irank << ", irank_head: " << mpi::irank_head << ", nprocs_head: " << mpi::nprocs_head << std::endl;

//  // message-passing test for a two-block grid
//  int test = mpi::irank * 100;
//  if (mpi::irank_block == 0) {
//
//    if (mpi::irank_head > 0)
//      MPI_Send(&test, 1, MPI_INT, 0, 1, mpi::comm_head);
//
//    else
//      MPI_Recv(&test, 1, MPI_INT, 1, 1, mpi::comm_head, mpi::status);
//
//  } // mpi::irank_block
//  std::cout << "Rank: " << std::setw(4) << mpi::irank << ", test value: " << std::setw(8) << test << std::endl;

  DEALLOCATE_1DPTR(ranks);

  mpi::wait_allothers();

  return;

} // group_cores_head



void label_cores_within_block(StructuredGrid *mygrid, StructuredBlock *myblock) {

  // within a block, label cores so that they can be easily accessed later

  int nshifts_core[DIM_MAX];

  // compute the neighbors' local ranks WITHOUT considering the block's decomposition
  nshifts_core[XI] = 1;
  nshifts_core[ETA] = myblock->num_cores_dir[XI];
  nshifts_core[ZETA] = myblock->num_cores_dir[XI] * myblock->num_cores_dir[ETA];
  //
  for (int idir = XI; idir < DIM_MAX; idir++) {

    mygrid->irank_next[idir][LEFT] = mpi::irank_block - nshifts_core[idir];
    mygrid->irank_next[idir][RIGHT] = mpi::irank_block + nshifts_core[idir];

  } // idir

  // adjust the neighbors' local ranks based upon the block's decomposition
  nshifts_core[XI] = myblock->num_cores_dir[XI] - 1;
  nshifts_core[ETA] = myblock->num_cores_dir[XI] * (myblock->num_cores_dir[ETA] - 1);
  nshifts_core[ZETA] = myblock->num_cores_dir[XI] * myblock->num_cores_dir[ETA] * (myblock->num_cores_dir[ZETA] - 1);
  //
  for (int idir = XI; idir < DIM_MAX; idir++) {

    // leftmost core
    if (mygrid->is_in_parent[idir] == myblock->is[idir]) {
      if (myblock->periodic[idir] == PERIODIC_PLANE)
        mygrid->irank_next[idir][LEFT] = mpi::irank_block + nshifts_core[idir];

      else
        mygrid->irank_next[idir][LEFT] = NONE; // now subject to either boundary condition or overset grid interpolation

    } // mygrid->is_in_parent[idir]

    // rightmost core
    if (mygrid->ie_in_parent[idir] == myblock->ie[idir]) {
      if (myblock->periodic[idir] == PERIODIC_PLANE)
        mygrid->irank_next[idir][RIGHT] = mpi::irank_block - nshifts_core[idir];

      else
        mygrid->irank_next[idir][RIGHT] = NONE; // now subject to either boundary condition or overset grid interpolation

    } // mygrid->je_in_parent[idir]
  } // idir

  //{
  //  std::cout << "Rank: " << mpi::irank << ", block: " << mygrid->id_parent << ", local rank: " << mpi::irank_block << ", next ranks: " 
  //            << mygrid->irank_next[XI][LEFT] << " " << mygrid->irank_next[XI][RIGHT] << ", "
  //            << mygrid->irank_next[ETA][LEFT] << " " << mygrid->irank_next[ETA][RIGHT] << ", "
  //            << mygrid->irank_next[ZETA][LEFT] << " " << mygrid->irank_next[ZETA][RIGHT] << std::endl;
  //}

  return;

} // label_cores_within_block



void StructuredGrid::override_local_indices(UserInput *myinput) {

  override_local_indices_4derivative(myinput);

  if (myinput->do_filter == TRUE)
    override_local_indices_4filter(myinput);

  return;

} // StructuredGrid::override_local_indices



void StructuredGrid::override_local_indices_4derivative(UserInput *myinput) {

  // by default, no adjustment on local indices
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is_4derivative[idir] = this->is[idir];
    this->ie_4derivative[idir] = this->ie[idir];

    this->iso_4derivative[idir] = this->iso[idir];
    this->ieo_4derivative[idir] = this->ieo[idir];

  } // idir

  // if a boundary does not have internal neighbors, adjust my indices for boundary derivatives
  for (int idir = XI; idir < DIM_MAX; idir++) {

    if (this->irank_next[idir][LEFT] == NONE && this->num_cells_dir[idir] > 1)
      this->is_4derivative[idir] += myinput->num_cells_ghost_finiteDifference;

    if (this->irank_next[idir][RIGHT] == NONE && this->num_cells_dir[idir] > 1)
      this->ie_4derivative[idir] -= myinput->num_cells_ghost_finiteDifference;

  } // idir

  for (int idir = XI; idir < myinput->num_dim; idir++) {

    this->iso_4derivative[idir] = this->is_4derivative[idir] - myinput->num_cells_ghost_finiteDifference;
    this->ieo_4derivative[idir] = this->ie_4derivative[idir] + myinput->num_cells_ghost_finiteDifference;

  } // idir

  return;

} // StructuredGrid::override_local_indices_4derivative



void StructuredGrid::override_local_indices_4filter(UserInput *myinput) {

  // by default, no adjustment on local indices
  for (int idir = XI; idir < DIM_MAX; idir++) {

    this->is_4filter[idir] = this->is[idir];
    this->ie_4filter[idir] = this->ie[idir];

    this->iso_4filter[idir] = this->iso[idir];
    this->ieo_4filter[idir] = this->ieo[idir];

  } // idir

  // if a boundary does not have internal neighbors, adjust my indices for boundary filters
  for (int idir = XI; idir < DIM_MAX; idir++) {

    if (this->irank_next[idir][LEFT] == NONE && this->num_cells_dir[idir] > 1)
      this->is_4filter[idir] += myinput->num_cells_ghost_filter;

    if (this->irank_next[idir][RIGHT] == NONE && this->num_cells_dir[idir] > 1)
      this->ie_4filter[idir] -= myinput->num_cells_ghost_filter;

  } // idir

  for (int idir = XI; idir < myinput->num_dim; idir++) {

    this->iso_4filter[idir] = this->is_4filter[idir] - myinput->num_cells_ghost_filter;
    this->ieo_4filter[idir] = this->ie_4filter[idir] + myinput->num_cells_ghost_filter;

  } // idir

  return;

} // StructuredGrid::override_local_indices_4filter



void additionalInit_boundary(UserInput *myinput, StructuredGrid *mygrid, StructuredBlock *block) {

  // do additional initialization or set-up for boundary treatment including boundary condition and buffer zone

  additionalInit_boundaryConditions(myinput, mygrid, block);
  additionalInit_bufferZones(myinput, mygrid, block);

  return;

} // additionalInit_boundary



void additionalInit_boundaryConditions(UserInput *myinput, StructuredGrid *mygrid, StructuredBlock *block) {



  return;

} // additionalInit_boundaryConditions



void additionalInit_bufferZones(UserInput *myinput, StructuredGrid *mygrid, StructuredBlock *block) {

  // check if the current block has any buffer zones within its grids; if not, return
  int num_bufferZones = mygrid->num_bufferZones;
  int num_bufferZones_sum = 0;
  MPI_Allreduce(&num_bufferZones, &num_bufferZones_sum, 1, MPI_INT, MPI_SUM, mpi::comm_block);
  //
//  std::cout << "Rank: " << mpi::irank << ", block: " << mygrid->id_parent << ", local rank: " << mpi::irank_block
//            << ", num of buffer zones: " << mygrid->num_bufferZones
//            << ", sum of num of buffer zones: " << num_bufferZones_sum
//            << std::endl;
  if (num_bufferZones_sum == 0)
    return;

  // the current block has at least one buffer zone
  int cur_block = mygrid->id_parent;
  int num_data_points = block[cur_block].nXi * block[cur_block].nEta * block[cur_block].nZeta;

  // loop over all buffer zones I have and get the xyz locations of buffer-zone boundaries in the block level
  for (int ibuffer = FIRST; ibuffer < num_bufferZones; ibuffer++) {

    StructuredBufferZone *bufferZone_cur = &(mygrid->bufferZone[ibuffer]);
    int idir_buf = bufferZone_cur->which_dir; // XI? ETA? ZETA?
    int iend_buf = bufferZone_cur->which_end; // LEFT? RIGHT?

    // initialize xyz locations on inflowing and outflowing sides of the parent buffer zone of mine
    for (int iside = IN; iside <= OUT; iside++) {

      StructuredBoundary *boundary_cur = &(bufferZone_cur->boundary[iside]);

      for (int k = boundary_cur->is[ZETA]; k <= boundary_cur->ie[ZETA]; k++) {
        int k_in_block = k - boundary_cur->is[ZETA] + boundary_cur->is_in_parent[ZETA];

        for (int j = boundary_cur->is[ETA]; j <= boundary_cur->ie[ETA]; j++) {
          int j_in_block = j - boundary_cur->is[ETA] + boundary_cur->is_in_parent[ETA];

          for (int i = boundary_cur->is[XI]; i <= boundary_cur->ie[XI]; i++) {
            int i_in_block = i - boundary_cur->is[XI] + boundary_cur->is_in_parent[XI];

            int lb = boundary_cur->idx1D(i, j, k); // 1-D index for boundary_cur
            int l0 = block[cur_block].idx1D(i_in_block, j_in_block, k_in_block); // 1-D index for xyz_block

            for (int idir = XDIR; idir < DIM_MAX; idir++)
              (boundary_cur->xyz[idir])[lb] = xyz_block[l0 + num_data_points * idir];

          } // i
        } // j
      } // k
    } // iside
  } // ibuffer

  // each buffer zone calculates the relative distance to the buffer zone boundaries
  double xyz_buf[DIM_MAX];
  double distance_2boundary[2]; // distances to inflow and outflow boundaries of the buffer zone
  for (int ibuffer = FIRST; ibuffer < num_bufferZones; ibuffer++) {

    StructuredBufferZone *bufferZone_cur = &(mygrid->bufferZone[ibuffer]);
    int idir_buf = bufferZone_cur->which_dir; // XI? ETA? ZETA?
    int iend_buf = bufferZone_cur->which_end; // LEFT? RIGHT?

    for (int k = bufferZone_cur->is[ZETA]; k <= bufferZone_cur->ie[ZETA]; k++) {
      int k_in_grid = k - bufferZone_cur->is[ZETA] + bufferZone_cur->is_in_parent[ZETA];

      for (int j = bufferZone_cur->is[ETA]; j <= bufferZone_cur->ie[ETA]; j++) {
        int j_in_grid = j - bufferZone_cur->is[ETA] + bufferZone_cur->is_in_parent[ETA];

        for (int i = bufferZone_cur->is[XI]; i <= bufferZone_cur->ie[XI]; i++) {
          int i_in_grid = i - bufferZone_cur->is[XI] + bufferZone_cur->is_in_parent[XI];

          int lb = bufferZone_cur->idx1D(i, j, k);
          int l0 = mygrid->idx1D(i_in_grid, j_in_grid, k_in_grid);

          // xyz locations of a point within the current buffer zone
          for (int idir = XDIR; idir < DIM_MAX; idir++)
            xyz_buf[idir] = mygrid->cell[l0].xyz[idir];

          for (int iside = IN; iside <= OUT; iside++) {

            // fetch a sponge boundary
            StructuredBoundary *boundary_cur = &(bufferZone_cur->boundary[iside]);
            distance_2boundary[iside] = 0.0;

            double distance_min = abs(DUMMY_DOUBLE);
            for (int lsb = 0; lsb < boundary_cur->num_cells; lsb++) {

              double distance = 0.0;
              for (int idir = XDIR; idir < DIM_MAX; idir++)
                distance += pow(xyz_buf[idir] - (boundary_cur->xyz[idir])[lsb], 2);

              distance_min = std::min(distance_min, distance);

            } // lsb

            distance_2boundary[iside] = sqrt(distance_min);

          } // iside

          // precompute the sponge strength, sigma * (relative-distance)^n
          double relative_distance = distance_2boundary[IN] / (distance_2boundary[IN] + distance_2boundary[OUT]);
          (bufferZone_cur->buffer_strength)[lb] = bufferZone_cur->buffer_constant * pow(relative_distance, bufferZone_cur->buffer_polynomial_order);

        } // i
      } // j
    } // k

//    for (int lb = 0; lb < bufferZone_cur->num_cells; lb++)
//      std::cout << "Block: " << cur_block << ", rank: " << mpi::irank
//                << ", ibuffer: " << ibuffer << ", idir: " << idir_buf << ", iend: " << iend_buf
//                << ", buffer strength: " << (bufferZone_cur->buffer_strength)[lb] << std::endl;;

  } // ibuffer

  return;

} // additionalInit_bufferZones

} // Geometry
