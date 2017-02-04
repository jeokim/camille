#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "plot3d.h"

namespace plot3d {

void pack_blockindex_2send(int *buf, Geometry::StructuredGrid *mygrid) {

  // order in which data are stored: js_XI, je_XI, js_ETA, je_ETA, js_ZETA, and je_ZETA

  for (int idir = XI; idir < DIM_MAX; idir++) {

    buf[    2 * idir] = mygrid->is_in_parent[idir]; // 2 since we have a pair of indices
    buf[1 + 2 * idir] = mygrid->ie_in_parent[idir];

  } // idir

  return;

} // pack_blockindex_2send



void pack_gridpoint_2send(double *dbuf, Geometry::StructuredGrid *mygrid) {

  // order in which data are stored: x(I,J,K = 0, ... ), y(I,J,K = 0, ... ), and z(I,J,K = 0, ... )

  for (int idir = XI; idir < DIM_MAX; idir++) {

    int num_points = 0; // counter for those who go into a 1-D array dbuf
    for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
      for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
        for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++) {

          int l0 = mygrid->idx1D(i, j, k);

          dbuf[mygrid->num_cells * idir + num_points] = mygrid->cell[l0].xyz[idir];
          num_points++;

        } // i
  } // idir

  return;

} // pack_gridpoint_2send



void pack_iblank_2send(int *ibuf, Geometry::StructuredGrid *mygrid) {

  int num_points = 0; // counter for those who go into a 1-D array ibuf
  for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
    for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
      for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++) {

        int l0 = mygrid->idx1D(i, j, k);

        ibuf[num_points] = mygrid->cell[l0].iblank;
        num_points++;

      } // i

  return;

} // pack_iblank_2send



void pack_solution_2send(double *dbuf, Geometry::StructuredGrid *mygrid, State *mystate) {

  // order in which data are stored: sol_0(I,J,K = 0, ... ), sol_1(I,J,K = 0, ... ), sol_2(I,J,K = 0, ... ), ...

  int num_vars = mystate->num_vars_sol;

  for (int ivar = 0; ivar < num_vars; ivar++) {

    int num_points = 0; // counter for those who go into a 1-D array dbuf
    for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
      for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
        for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++) {

          int l0 = mygrid->idx1D(i, j, k);

          dbuf[mygrid->num_cells * ivar + num_points] = (mystate->sol[ivar])[l0];
          num_points++;

        } // i
  } // ivar

  return;

} // pack_solution_2send



void pack_function_2send(double *dbuf, Geometry::StructuredGrid *mygrid, int num_vars, double **func) {

  // order in which data are stored: func_0(I,J,K = 0, ... ), func_1(I,J,K = 0, ... ), func_2(I,J,K = 0, ... ), ...

  for (int ivar = 0; ivar < num_vars; ivar++) {

    int num_points = 0; // counter for those who go into a 1-D array dbuf
    for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
      for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
        for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++) {

          int l0 = mygrid->idx1D(i, j, k);

          dbuf[mygrid->num_cells * ivar + num_points] = (func[ivar])[l0];
          num_points++;

        } // i
  } // ivar

  return;

} // pack_function_2send



void send_blockindex_2headrank(Geometry::StructuredGrid *mygrid, int *idxs_in_parent, int *idxe_in_parent) {

  int num_grids;
  //
  int count;
  int *buf;

  num_grids = mygrid->num_ourkinds;

  count = 2 * DIM_MAX; //  2 because there is a pair of is and ie (or js and je)
  buf = new int[count];

  // message passing within a block to deliver index information to a head rank
  for (int id_local = 0; id_local < num_grids; id_local++) {

    for (int i = 0; i < count; i++)
      buf[i] = -1; // some dummy number

    if (mygrid->id_local == id_local) // every core packs its own indices
      pack_blockindex_2send(buf, mygrid);

    if (id_local > 0) {

      if (mygrid->id_local == id_local) // non-head ranks send
        MPI_Send(buf, count, MPI_INT, 0, id_local, mpi::comm_block);

      else if (mygrid->id_local == 0) // a head rank always receives
        MPI_Recv(buf, count, MPI_INT, id_local, id_local, mpi::comm_block, mpi::status);

    } // id_local

    if (mpi::irank_block == 0) { // a head rank stores
      for (int idir = XI; idir < DIM_MAX; idir++) {

        idxs_in_parent[id_local * DIM_MAX + idir] = buf[    2 * idir];
        idxe_in_parent[id_local * DIM_MAX + idir] = buf[1 + 2 * idir];

      } // idir
    } // mpi::irank_block
  } // id_local

  // clean up
  DEALLOCATE_1DPTR(buf);
  mpi::wait_allothers();
  //
  //if (mpi::irank_block == 0)
  //    for (int id_local = 0; id_local < num_grids; id_local++)
  //      std::cout << "Rank: " << mpi::irank <<", id_local: " << id_local << ", is-ie (parent): " 
  //                                               << idxs_in_parent[id_local * DIM_MAX + XI] << " " << idxe_in_parent[id_local * DIM_MAX + XI] << " "
  //                                               << idxs_in_parent[id_local * DIM_MAX + ETA] << " " << idxe_in_parent[id_local * DIM_MAX + ETA] << " "
  //                                               << idxs_in_parent[id_local * DIM_MAX + ZETA] << " " << idxe_in_parent[id_local * DIM_MAX + ZETA] << " "
  //                                               << std::endl;

  return;

} // send_blockindex_2headrank



void read_grid_header(std::string filename, int &num_blocks_in_file, int *&num_cells_in) {

  if (mpi::irank == 0) {

    std::ifstream ifs;

    ifs.open(cstr_to_constchar(filename), std::ifstream::binary | std::ifstream::in);
    if (!ifs.is_open())
      mpi::graceful_exit("The PLOT3D grid file does not exist.");

    ifs.read(reinterpret_cast<char *>(&num_blocks_in_file), sizeof num_blocks_in_file);
    ALLOCATE1D_INT_2ARG(num_cells_in, num_blocks_in_file, DIM_MAX);

    int count = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
      for (int idir = XI; idir < DIM_MAX; idir++) {

        ifs.read(reinterpret_cast<char *>(&num_cells_in[count]), sizeof num_cells_in[count]);
        count++;

      } // idir

    ifs.close();

  } // mpi::irank

  mpi::wait_allothers();

  // now broadcast
  MPI_Bcast(&num_blocks_in_file, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (mpi::irank != 0)
    ALLOCATE1D_INT_2ARG(num_cells_in, num_blocks_in_file, DIM_MAX);

  MPI_Bcast(num_cells_in, num_blocks_in_file * DIM_MAX, MPI_INT, 0, MPI_COMM_WORLD);

  return;

} // read_grid_header



void read_grid_header_fileKeptOpen(std::string filename, int &num_blocks_in_file, int *&num_cells_in, std::ifstream &ifs) {

  // this function does basically the same thing as read_grid_header
  // one difference is that file is opened and not closed for future use

  if (mpi::irank == 0) {

    //std::ifstream ifs;

    ifs.open(cstr_to_constchar(filename), std::ifstream::binary | std::ifstream::in);
    if (!ifs.is_open())
      mpi::graceful_exit("The PLOT3D grid file does not exist.");

    ifs.read(reinterpret_cast<char *>(&num_blocks_in_file), sizeof num_blocks_in_file);
    ALLOCATE1D_INT_2ARG(num_cells_in, num_blocks_in_file, DIM_MAX);

    int count = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
      for (int idir = XI; idir < DIM_MAX; idir++) {

        ifs.read(reinterpret_cast<char *>(&num_cells_in[count]), sizeof num_cells_in[count]);
        count++;

      } // idir

    //ifs.close();

  } // mpi::irank

  mpi::wait_allothers();

  // now broadcast
  MPI_Bcast(&num_blocks_in_file, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (mpi::irank != 0)
    ALLOCATE1D_INT_2ARG(num_cells_in, num_blocks_in_file, DIM_MAX);

  MPI_Bcast(num_cells_in, num_blocks_in_file * DIM_MAX, MPI_INT, 0, MPI_COMM_WORLD);

  return;

} // read_grid_header_fileKeptOpen



void read_solution_header(std::string filename, int &num_blocks_in_file, int *&num_cells_in) {

  // PLOT3D solution q file has the same header as PLOT3D grid xyz file
  read_grid_header(filename, num_blocks_in_file, num_cells_in);

  return;

} // read_solution_header



void read_solution_header_fileKeptOpen(std::string filename, int &num_blocks_in_file, int *&num_cells_in, std::ifstream &ifs) {

  // PLOT3D solution q file has the same header as PLOT3D grid xyz file
  read_grid_header_fileKeptOpen(filename, num_blocks_in_file, num_cells_in, ifs);

  return;

} // read_solution_header_fileKeptOpen



void stdout_grid_size(std::string filename, int whichPLOT3Dfile) {

  int num_blocks_in_file;
  int *num_cells_dir;
  std::string str_stdout;
  std::stringstream str_counter;

  switch ( whichPLOT3Dfile ) {
  case PLOT3D_GRID:

    str_stdout = "The grid file contains ";
    read_grid_header(filename, num_blocks_in_file, num_cells_dir);

    break;

  case PLOT3D_SOLUTION:

    str_stdout = "The solution file contains ";
    read_solution_header(filename, num_blocks_in_file, num_cells_dir);

    break;

  default:

    mpi::graceful_exit("Unknown type of PLOT3D file: either grid, solution, or function.");

    break;

  } // whichPLOT3Dfile

  str_counter << std::setw(6) << num_blocks_in_file;
  str_stdout += str_counter.str() + " blocks.";
  MESSAGE_STDOUT(str_stdout);

  MESSAGE_STDOUT("");
  MESSAGE_STDOUT("block-ID     points-in-XI    points-in-ETA   points-in-ZETA            total");
  MESSAGE_STDOUT("----------------------------------------------------------------------------");

  int count = 0;
  int num_cells_total = 0;
  for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

    int num_cells_per_block = num_cells_dir[count]*num_cells_dir[count+1]*num_cells_dir[count+2];
    num_cells_total += num_cells_per_block;

    str_counter.str("");
    str_counter << std::setw(8) << iblock;
    str_stdout = str_counter.str();

    for (int idim = XI; idim < DIM_MAX; idim++) {

      str_counter.str("");
      str_counter << std::setw(17) << num_cells_dir[count++];
      str_stdout += str_counter.str();

    } // idim

    str_counter.str("");
    str_counter << std::setw(17) << num_cells_per_block;
    str_stdout += str_counter.str();

    MESSAGE_STDOUT(str_stdout);

  } // iblock
  MESSAGE_STDOUT("----------------------------------------------------------------------------");

  str_counter.str("");
  str_counter << std::setw(20) << num_cells_total/1000000.0;
  str_stdout = "Total number of points is " + str_counter.str() + " million.";
  MESSAGE_STDOUT(str_stdout);
  MESSAGE_STDOUT("");

  num_blocks_in_file = 0;
  delete[] num_cells_dir;

  return;

} // stdout_grid_size



void write_grid_header(std::string filename, Geometry::StructuredBlock *block) {

  if (mpi::irank == 0) {

    int cur_block = 0; // a global head rank is always in block-0
    int num_blocks = block[cur_block].num_ourkinds; // every block has the same value at block[i].num_ourkinds
    std::ofstream ofs;

    ofs.open(cstr_to_constchar(filename), std::ofstream::binary | std::ofstream::trunc);
    ofs.write(reinterpret_cast<const char*>(&num_blocks), sizeof num_blocks);

    for (int iblock = 0; iblock < num_blocks; iblock++)
      for (int idir = XI; idir < DIM_MAX; idir++)
        ofs.write(reinterpret_cast<const char*>(&block[iblock].num_cells_dir[idir]), sizeof block[iblock].num_cells_dir[idir]);

    ofs.close();

  } // mpi::irank

  mpi::wait_allothers();

  return;

} // write_grid_header



void write_solution_header(std::string filename, Geometry::StructuredBlock *block) {

  // PLOT3D solution q file has the same header as PLOT3D grid xyz file
  write_grid_header(filename, block);

  return;

} // write_solution_header



void write_function_header(std::string filename, Geometry::StructuredBlock *block, int num_vars) {

  if (mpi::irank == 0) {

    int cur_block = 0; // a global head rank is always in block-0
    int num_blocks = block[cur_block].num_ourkinds; // every block has the same value at block[i].num_ourkinds
    std::ofstream ofs;

    ofs.open(cstr_to_constchar(filename), std::ofstream::binary | std::ofstream::trunc);
    ofs.write(reinterpret_cast<const char*>(&num_blocks), sizeof num_blocks);

    for (int iblock = 0; iblock < num_blocks; iblock++) {

      for (int idir = XI; idir < DIM_MAX; idir++)
        ofs.write(reinterpret_cast<const char*>(&block[iblock].num_cells_dir[idir]), sizeof block[iblock].num_cells_dir[idir]);

      ofs.write(reinterpret_cast<const char*>(&num_vars), sizeof num_vars); // this line is an only difference with write_grid_headern

    } // iblock

    ofs.close();

  } // mpi::irank

  mpi::wait_allothers();

  return;

} // write_function_header



void read_grid_serialIO(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  int cur_block, num_blocks, num_points_cur_block;
  double *data_block;
  int num_vars_2read;
  std::string filename;
  //
  int num_blocks_in_file;

  // useful variables
  cur_block = mygrid->id_parent; // a block ID to which this core (or equivalently, grid) belongs
  num_blocks = block[cur_block].num_ourkinds; // number of blocks is shared by every block
  num_vars_2read = DIM_MAX; // DIM_MAX since grid points are written
  filename = myinput->file_grid_in;

  // allocate
  num_points_cur_block = block[cur_block].num_cells;
  data_block = new double[num_points_cur_block * num_vars_2read];

  // read my block's grid data
  for (int iblock = 0; iblock < num_blocks; iblock++) { // read in the order of block ID
    if (cur_block == iblock) { // to avoid simultaneous excessive access by too many cores, read only if I belong to this block

      int ibuf;
      double dbuf;

      std::ifstream ifs;

      ifs.open(cstr_to_constchar(filename), std::ifstream::binary | std::ifstream::in);
      if (!ifs.is_open())
        mpi::graceful_exit("The PLOT3D file named " + filename + " does not exist.");

      // read and throw away the file header information
      ifs.read(reinterpret_cast<char *>(&num_blocks_in_file), sizeof num_blocks_in_file);
      for (int iblock_in_file = 0; iblock_in_file < num_blocks_in_file; iblock_in_file++)
        for (int idir = XI; idir < DIM_MAX; idir++)
          ifs.read(reinterpret_cast<char *>(&ibuf), sizeof ibuf);

      // read and throw away previous blocks' data
      for (int iblock_former = 0; iblock_former < cur_block; iblock_former++) {

        for (int ivar = 0; ivar < num_vars_2read; ivar++)
          for (int k = 0; k < block[iblock_former].nZeta; k++)
            for (int j = 0; j < block[iblock_former].nEta; j++)
              for (int i = 0; i < block[iblock_former].nXi; i++)
                ifs.read(reinterpret_cast<char *>(&dbuf), sizeof dbuf);

      } // iblock_former

      // now read my block's
      int num_data_points = block[iblock].nXi * block[iblock].nEta * block[iblock].nZeta;
      for (int ivar = 0; ivar < num_vars_2read; ivar++) {
        for (int k = 0; k < block[iblock].nZeta; k++)
          for (int j = 0; j < block[iblock].nEta; j++)
            for (int i = 0; i < block[iblock].nXi; i++) {

              int l0 = block[iblock].idx1D(i, j, k);

              ifs.read(reinterpret_cast<char *>(&data_block[l0 + num_data_points * ivar]), sizeof data_block[0]);

            } // i
      } // ivar
      ifs.close();

    } // cur_block
    mpi::wait_allothers();

  } // iblock

  // each core has grid-point data for the block to which it belongs; now only take my portion
  // note that ghost cells are assigned their coordinate locations as well
  int num_data_points = block[cur_block].nXi * block[cur_block].nEta * block[cur_block].nZeta;
  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_parent = mygrid->iso_in_parent[ZETA] + k;
    k_in_parent = std::max(k_in_parent, 0); // ghost cell locations existing outside of block are clipped
    k_in_parent = std::min(k_in_parent, block[cur_block].nZeta - 1); // ghost cell locations existing outside of block are clipped

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_parent = mygrid->iso_in_parent[ETA] + j;
      j_in_parent = std::max(j_in_parent, 0); // ghost cell locations existing outside of block are clipped
      j_in_parent = std::min(j_in_parent, block[cur_block].nEta - 1); // ghost cell locations existing outside of block are clipped

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_parent = mygrid->iso_in_parent[XI] + i;
        i_in_parent = std::max(i_in_parent, 0); // ghost cell locations existing outside of block are clipped
        i_in_parent = std::min(i_in_parent, block[cur_block].nXi - 1); // ghost cell locations existing outside of block are clipped

        int l0 = mygrid->idx1D(i, j, k);
        int lb = block[cur_block].idx1D(i_in_parent, j_in_parent, k_in_parent);

        for (int idir = XDIR; idir < DIM_MAX; idir++)
          mygrid->cell[l0].xyz[idir] = data_block[lb + num_data_points * idir];

      } // i
    } // j
  } // k

  // keep my block's xyz data
  Geometry::fill_in_xyz_block(data_block, num_points_cur_block * num_vars_2read);

  DEALLOCATE_1DPTR(data_block);

//  // extrapolate non-internal ghost cells for cores near block boundaries
//  // in XI
//  if (mygrid->irank_next[XI][LEFT] == NONE && mygrid->num_cells_dir[XI] > 1) {
//    for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
//      for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
//        for (int i_stencil = 1; i_stencil < mygrid->num_cells_ghost; i_stencil++) {
//
//          int i = mygrid->is[XI] - i_stencil;
//          int l0 = mygrid->idx1D(i, j, k);
//
//          i = mygrid->is[XI];
//          int l1 = mygrid->idx1D(i, j, k);
//
//          i = mygrid->is[XI] + i_stencil;
//          int l2 = mygrid->idx1D(i, j, k);
//
//          for (int idir = XDIR; idir < DIM_MAX; idir++)
//            mygrid->cell[l0].xyz[idir] = 2.0 * mygrid->cell[l1].xyz[idir] - mygrid->cell[l2].xyz[idir];
//
//        } // i_stencil
//  } // mygrid->irank_next[XI][LEFT]
//  //
//  if (mygrid->irank_next[XI][RIGHT] == NONE && mygrid->num_cells_dir[XI] > 1) {
//    for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
//      for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
//        for (int i_stencil = 1; i_stencil < mygrid->num_cells_ghost; i_stencil++) {
//
//          int i = mygrid->is[XI] + i_stencil;
//          int l0 = mygrid->idx1D(i, j, k);
//
//          i = mygrid->is[XI];
//          int l1 = mygrid->idx1D(i, j, k);
//
//          i = mygrid->is[XI] - i_stencil;
//          int l2 = mygrid->idx1D(i, j, k);
//
//          for (int idir = XDIR; idir < DIM_MAX; idir++)
//            mygrid->cell[l0].xyz[idir] = 2.0 * mygrid->cell[l1].xyz[idir] - mygrid->cell[l2].xyz[idir];
//
//        } // i_stencil
//  } // mygrid->irank_next[XI][RIGHT]
//
//  // in ETA
//  if (mygrid->irank_next[ETA][LEFT] == NONE && mygrid->num_cells_dir[ETA] > 1) {
//    for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
//      for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++)
//        for (int j_stencil = 1; j_stencil < mygrid->num_cells_ghost; j_stencil++) {
//
//          int j = mygrid->is[ETA] - j_stencil;
//          int l0 = mygrid->idx1D(i, j, k);
//
//          j = mygrid->is[ETA];
//          int l1 = mygrid->idx1D(i, j, k);
//
//          j = mygrid->is[ETA] + j_stencil;
//          int l2 = mygrid->idx1D(i, j, k);
//
//          for (int idir = XDIR; idir < DIM_MAX; idir++)
//            mygrid->cell[l0].xyz[idir] = 2.0 * mygrid->cell[l1].xyz[idir] - mygrid->cell[l2].xyz[idir];
//
//        } // j_stencil
//  } // mygrid->irank_next[ETA][LEFT]
//  //
//  if (mygrid->irank_next[ETA][RIGHT] == NONE && mygrid->num_cells_dir[ETA] > 1) {
//    for (int k = mygrid->is[ZETA]; k < mygrid->ie[ZETA] + 1; k++)
//      for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++)
//        for (int j_stencil = 1; j_stencil < mygrid->num_cells_ghost; j_stencil++) {
//
//          int j = mygrid->is[ETA] + j_stencil;
//          int l0 = mygrid->idx1D(i, j, k);
//
//          j = mygrid->is[ETA];
//          int l1 = mygrid->idx1D(i, j, k);
//
//          j = mygrid->is[ETA] - j_stencil;
//          int l2 = mygrid->idx1D(i, j, k);
//
//          for (int idir = XDIR; idir < DIM_MAX; idir++)
//            mygrid->cell[l0].xyz[idir] = 2.0 * mygrid->cell[l1].xyz[idir] - mygrid->cell[l2].xyz[idir];
//
//        } // j_stencil
//  } // mygrid->irank_next[ETA][RIGHT]
//
//  // in ZETA
//  if (mygrid->irank_next[ZETA][LEFT] == NONE && mygrid->num_cells_dir[ZETA] > 1) {
//    for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
//      for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++)
//        for (int k_stencil = 1; k_stencil < mygrid->num_cells_ghost; k_stencil++) {
//
//          int k = mygrid->is[ZETA] - k_stencil;
//          int l0 = mygrid->idx1D(i, j, k);
//
//          k = mygrid->is[ZETA];
//          int l1 = mygrid->idx1D(i, j, k);
//
//          k = mygrid->is[ZETA] + k_stencil;
//          int l2 = mygrid->idx1D(i, j, k);
//
//          for (int idir = XDIR; idir < DIM_MAX; idir++)
//            mygrid->cell[l0].xyz[idir] = 2.0 * mygrid->cell[l1].xyz[idir] - mygrid->cell[l2].xyz[idir];
//
//        } // k_stencil
//  } // mygrid->irank_next[ZETA][LEFT]
//  //
//  if (mygrid->irank_next[ZETA][RIGHT] == NONE && mygrid->num_cells_dir[ZETA] > 1) {
//    for (int j = mygrid->is[ETA]; j < mygrid->ie[ETA] + 1; j++)
//      for (int i = mygrid->is[XI]; i < mygrid->ie[XI] + 1; i++)
//        for (int k_stencil = 1; k_stencil < mygrid->num_cells_ghost; k_stencil++) {
//
//          int k = mygrid->is[ZETA] + k_stencil;
//          int l0 = mygrid->idx1D(i, j, k);
//
//          k = mygrid->is[ZETA];
//          int l1 = mygrid->idx1D(i, j, k);
//
//          k = mygrid->is[ZETA] - k_stencil;
//          int l2 = mygrid->idx1D(i, j, k);
//
//          for (int idir = XDIR; idir < DIM_MAX; idir++)
//            mygrid->cell[l0].xyz[idir] = 2.0 * mygrid->cell[l1].xyz[idir] - mygrid->cell[l2].xyz[idir];
//
//        } // k_stencil
//  } // mygrid->irank_next[ZETA][RIGHT]

  return;

} // read_grid_serialIO



void minmax_in_gridFile(int num_blocks_in_grid, int *num_cells_in, int num_vars, double *&minmax, std::ifstream &ifs) {

  // this function reads a PLOT3D grid file after the header (i.e. ifs is assumed to have read the header)
  // min and max of each variable are obtained

  if (!ifs.is_open())
    mpi::graceful_exit("The requested PLOT3D file is not opened.");
  if (num_vars != DIM_MAX)
    mpi::graceful_exit("PLOT3D grid file should have 3 variables (x, y, and z).");

  for (int ivar = 0; ivar < num_vars; ivar++) {

    minmax[2*ivar+MINIMUM] =  DUMMY_LARGE;
    minmax[2*ivar+MAXIMUM] = -DUMMY_LARGE;

  } // ivar

  int count = 0;
  for (int iblock = 0; iblock < num_blocks_in_grid; iblock++) {

    double dbuf;
    for (int ivar = 0; ivar < num_vars; ivar++)
      for (int k = 0; k < num_cells_in[count+2]; k++)
        for (int j = 0; j < num_cells_in[count+1]; j++)
          for (int i = 0; i < num_cells_in[count]; i++) {

            ifs.read(reinterpret_cast<char *>(&dbuf), sizeof dbuf);
            minmax[2*ivar+MINIMUM] = std::min(minmax[2*ivar+MINIMUM], dbuf);
            minmax[2*ivar+MAXIMUM] = std::max(minmax[2*ivar+MAXIMUM], dbuf);

          } // i
    count += DIM_MAX;
  } // iblock

  return;

} // minmax_in_gridFile



void read_solution_serialIO(std::string filename, UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate) {

  int cur_block, num_blocks, num_points_cur_block;
  double *data_block;
  int num_vars_2read;
  //
  int num_blocks_in_file;
  double tau[4];

  // useful variables
  cur_block = mygrid->id_parent; // a block ID to which this core (or equivalently, grid) belongs
  num_blocks = block[cur_block].num_ourkinds; // number of blocks is shared by every block
  num_vars_2read = mystate->num_vars_sol;

  // allocate
  num_points_cur_block = block[cur_block].num_cells;
  data_block = new double[num_points_cur_block * num_vars_2read];

  // read my block's grid data
  for (int iblock = 0; iblock < num_blocks; iblock++) { // read in the order of block ID
    if (cur_block == iblock) { // to avoid simultaneous excessive access by too many cores, read only if I belong to this block

      int ibuf[DIM_MAX];
      double dbuf;

      std::ifstream ifs;

      ifs.open(cstr_to_constchar(filename), std::ifstream::binary | std::ifstream::in);
      if (!ifs.is_open())
        mpi::graceful_exit("The PLOT3D file named " + filename + " does not exist.");

      // read and throw away the file header information
      ifs.read(reinterpret_cast<char *>(&num_blocks_in_file), sizeof num_blocks_in_file);
      for (int iblock_in_file = 0; iblock_in_file < num_blocks_in_file; iblock_in_file++) {

        for (int idir = XI; idir < DIM_MAX; idir++)
          ifs.read(reinterpret_cast<char *>(&ibuf[idir]), sizeof ibuf[idir]);

        // check if the solution file contains a consistent number of data as the corresponding grid file
        assert( ibuf[XI] == block[iblock_in_file].nXi );
        assert( ibuf[ETA] == block[iblock_in_file].nEta );
        assert( ibuf[ZETA] == block[iblock_in_file].nZeta );

      } // iblock_in_file

      // read and throw away previous blocks' data
      for (int iblock_former = 0; iblock_former < cur_block; iblock_former++) {

        // specific to PLOT3D solution file
        ifs.read(reinterpret_cast<char *>(&tau), sizeof tau);

        for (int ivar = 0; ivar < num_vars_2read; ivar++)
          for (int k = 0; k < block[iblock_former].nZeta; k++)
            for (int j = 0; j < block[iblock_former].nEta; j++)
              for (int i = 0; i < block[iblock_former].nXi; i++)
                ifs.read(reinterpret_cast<char *>(&dbuf), sizeof dbuf);

      } // iblock_former

      // now read my block's
      // specific to PLOT3D solution file
      ifs.read(reinterpret_cast<char *>(&tau), sizeof tau);
      mystate->time_step = static_cast<int>(tau[2]); // get the previous run's time step
      mystate->time_step_lastrun = mystate->time_step;
      mystate->time_sol = tau[3]; // get the previous run's solution time
      mystate->time_sol_lastrun = mystate->time_sol;
      //
      int num_data_points = block[iblock].nXi * block[iblock].nEta * block[iblock].nZeta;
      for (int ivar = 0; ivar < num_vars_2read; ivar++) {
        for (int k = 0; k < block[iblock].nZeta; k++)
          for (int j = 0; j < block[iblock].nEta; j++)
            for (int i = 0; i < block[iblock].nXi; i++) {

              int l0 = block[iblock].idx1D(i, j, k);

              ifs.read(reinterpret_cast<char *>(&data_block[l0 + num_data_points * ivar]), sizeof data_block[0]);

            } // i
      } // ivar
      ifs.close();

    } // cur_block
    mpi::wait_allothers();

  } // iblock

  // each core has grid-point data for the block to which it belongs; now only take my portion
  // note that ghost cells are assigned their solution variables as well
  int num_data_points = block[cur_block].nXi * block[cur_block].nEta * block[cur_block].nZeta;
  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_parent = mygrid->iso_in_parent[ZETA] + k;
    k_in_parent = std::max(k_in_parent, 0); // ghost cell locations existing outside of block are clipped
    k_in_parent = std::min(k_in_parent, block[cur_block].nZeta - 1); // ghost cell locations existing outside of block are clipped

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_parent = mygrid->iso_in_parent[ETA] + j;
      j_in_parent = std::max(j_in_parent, 0); // ghost cell locations existing outside of block are clipped
      j_in_parent = std::min(j_in_parent, block[cur_block].nEta - 1); // ghost cell locations existing outside of block are clipped

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_parent = mygrid->iso_in_parent[XI] + i;
        i_in_parent = std::max(i_in_parent, 0); // ghost cell locations existing outside of block are clipped
        i_in_parent = std::min(i_in_parent, block[cur_block].nXi - 1); // ghost cell locations existing outside of block are clipped

        int l0 = mygrid->idx1D(i, j, k);
        int lb = block[cur_block].idx1D(i_in_parent, j_in_parent, k_in_parent);

        for (int ivar = 0; ivar < num_vars_2read; ivar++)
          (mystate->sol[ivar])[l0] = data_block[lb + num_data_points * ivar];

      } // i
    } // j
  } // k

  DEALLOCATE_1DPTR(data_block);

  return;

} // read_solution_serialIO



void minmax_in_solutionFile(int num_blocks_in_grid, int *num_cells_in, int num_vars, double *&tau_in, double *&minmax, std::ifstream &ifs) {

  // this function reads a PLOT3D solution file after the header (i.e. ifs is assumed to have read the header)
  // min and max of each variable are obtained

  if (!ifs.is_open())
    mpi::graceful_exit("The requested PLOT3D file is not opened.");
  if (num_vars != DIM_MAX+2)
    mpi::graceful_exit("PLOT3D solution file should have 5 solution variables.");

  for (int ivar = 0; ivar < num_vars; ivar++) {

    minmax[2*ivar+MINIMUM] =  DUMMY_LARGE;
    minmax[2*ivar+MAXIMUM] = -DUMMY_LARGE;

  } // ivar

  int count = 0;
  double tau[4];
  for (int iblock = 0; iblock < num_blocks_in_grid; iblock++) {
    ifs.read(reinterpret_cast<char *>(&tau), sizeof tau);
    for (int ivar = 0; ivar < 4; ivar++)
      tau_in[ivar] = tau[ivar];

    double dbuf;
    for (int ivar = 0; ivar < num_vars; ivar++)
      for (int k = 0; k < num_cells_in[count+2]; k++)
        for (int j = 0; j < num_cells_in[count+1]; j++)
          for (int i = 0; i < num_cells_in[count]; i++) {

            ifs.read(reinterpret_cast<char *>(&dbuf), sizeof dbuf);
            minmax[2*ivar+MINIMUM] = std::min(minmax[2*ivar+MINIMUM], dbuf);
            minmax[2*ivar+MAXIMUM] = std::max(minmax[2*ivar+MAXIMUM], dbuf);

          } // i
    count += DIM_MAX;
  } // iblock

  return;

} // minmax_in_solutionFile



void read_function_serialIO(std::string filename, UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, int num_vars, double **var) {

  int cur_block, num_blocks, num_points_cur_block;
  double *data_block;
  int num_vars_2read;
  //
  int num_blocks_in_file;

  // useful variables
  cur_block = mygrid->id_parent; // a block ID to which this core (or equivalently, grid) belongs
  num_blocks = block[cur_block].num_ourkinds; // number of blocks is shared by every block
  num_vars_2read = num_vars;

  // allocate
  num_points_cur_block = block[cur_block].num_cells;
  data_block = new double[num_points_cur_block * num_vars_2read];

  // read my block's grid data
  for (int iblock = 0; iblock < num_blocks; iblock++) { // read in the order of block ID
    if (cur_block == iblock) { // to avoid simultaneous excessive access by too many cores, read only if I belong to this block

      int ibuf[DIM_MAX + 1];
      double dbuf;

      std::ifstream ifs;

      ifs.open(cstr_to_constchar(filename), std::ifstream::binary | std::ifstream::in);
      if (!ifs.is_open())
        mpi::graceful_exit("The PLOT3D file named " + filename + " does not exist.");

      // read and throw away the file header information
      ifs.read(reinterpret_cast<char *>(&num_blocks_in_file), sizeof num_blocks_in_file);
      for (int iblock_in_file = 0; iblock_in_file < num_blocks_in_file; iblock_in_file++) {

        for (int idir = XI; idir < DIM_MAX + 1; idir++)
          ifs.read(reinterpret_cast<char *>(&ibuf[idir]), sizeof ibuf[idir]);

        // check if the solution file contains a consistent number of data as the corresponding grid file
        assert( ibuf[XI] == block[iblock_in_file].nXi );
        assert( ibuf[ETA] == block[iblock_in_file].nEta );
        assert( ibuf[ZETA] == block[iblock_in_file].nZeta );
        assert( ibuf[DIM_MAX] == num_vars );

      } // iblock_in_file

      // read and throw away previous blocks' data
      for (int iblock_former = 0; iblock_former < cur_block; iblock_former++) {

        for (int ivar = 0; ivar < num_vars_2read; ivar++)
          for (int k = 0; k < block[iblock_former].nZeta; k++)
            for (int j = 0; j < block[iblock_former].nEta; j++)
              for (int i = 0; i < block[iblock_former].nXi; i++)
                ifs.read(reinterpret_cast<char *>(&dbuf), sizeof dbuf);

      } // iblock_former

      // now read my block's
      int num_data_points = block[iblock].nXi * block[iblock].nEta * block[iblock].nZeta;
      for (int ivar = 0; ivar < num_vars_2read; ivar++) {
        for (int k = 0; k < block[iblock].nZeta; k++)
          for (int j = 0; j < block[iblock].nEta; j++)
            for (int i = 0; i < block[iblock].nXi; i++) {

              int l0 = block[iblock].idx1D(i, j, k);

              ifs.read(reinterpret_cast<char *>(&data_block[l0 + num_data_points * ivar]), sizeof data_block[0]);

            } // i
      } // ivar
      ifs.close();

    } // cur_block
    mpi::wait_allothers();

  } // iblock

  // each core has grid-point data for the block to which it belongs; now only take my portion
  // note that ghost cells are assigned their solution variables as well
  int num_data_points = block[cur_block].nXi * block[cur_block].nEta * block[cur_block].nZeta;
  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    int k_in_parent = mygrid->iso_in_parent[ZETA] + k;
    k_in_parent = std::max(k_in_parent, 0); // ghost cell locations existing outside of block are clipped
    k_in_parent = std::min(k_in_parent, block[cur_block].nZeta - 1); // ghost cell locations existing outside of block are clipped

    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      int j_in_parent = mygrid->iso_in_parent[ETA] + j;
      j_in_parent = std::max(j_in_parent, 0); // ghost cell locations existing outside of block are clipped
      j_in_parent = std::min(j_in_parent, block[cur_block].nEta - 1); // ghost cell locations existing outside of block are clipped

      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int i_in_parent = mygrid->iso_in_parent[XI] + i;
        i_in_parent = std::max(i_in_parent, 0); // ghost cell locations existing outside of block are clipped
        i_in_parent = std::min(i_in_parent, block[cur_block].nXi - 1); // ghost cell locations existing outside of block are clipped

        int l0 = mygrid->idx1D(i, j, k);
        int lb = block[cur_block].idx1D(i_in_parent, j_in_parent, k_in_parent);

        for (int ivar = 0; ivar < num_vars_2read; ivar++)
          (var[ivar])[l0] = data_block[lb + num_data_points * ivar];

      } // i
    } // j
  } // k

  DEALLOCATE_1DPTR(data_block);

  return;

} // read_function_serialIO



void write_grid_serialIO(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  int cur_block, num_blocks, num_points_cur_block;
  int num_grids;
  int num_vars_2write;
  int *idxs_in_parent, *idxe_in_parent;
  double *data_block;
  int *iblank_block;
  std::string filename;
  //
  int count;
  double *dbuf;
  int *ibuf;

  // useful variables
  cur_block = mygrid->id_parent; // a block ID to which this core (or equivalently, grid) belongs
  num_blocks = block[cur_block].num_ourkinds; // number of blocks is shared by every block
  num_grids = mygrid->num_ourkinds;
  num_vars_2write = DIM_MAX; // DIM_MAX since grid points are written
  filename = myinput->file_grid;

  // first, write a header for grid
  write_grid_header(filename, block);

  // second, send block indices to a head rank of each block
  // only a head rank in a block allocates to store
  if (mpi::irank_block == 0) {

    idxs_in_parent = new int[num_grids * DIM_MAX]; // a grid's starting index in the parent block
    idxe_in_parent = new int[num_grids * DIM_MAX]; // a grid's ending index in the parent block

  } // mpi::irank_block
  //
  send_blockindex_2headrank(mygrid, idxs_in_parent, idxe_in_parent);
  mpi::wait_allothers();

  // third, grid point information
  // array for local grid point data
  count = mygrid->num_cells * num_vars_2write;
  num_points_cur_block = block[cur_block].num_cells;
  dbuf = new double[count];
  ibuf = new int[mygrid->num_cells];
  if (mpi::irank_block == 0) { // only a head rank in a block allocates

    data_block = new double[num_points_cur_block * num_vars_2write];
    iblank_block = new int[num_points_cur_block];

  } // mpi::irank_block
  //
  // message passing within a block to deliver grid point information to a head rank
  for (int id_local = 0; id_local < num_grids; id_local++) {

    if (mygrid->id_local == id_local) {

      pack_gridpoint_2send(dbuf, mygrid); // every core packs grid point information
      pack_iblank_2send(ibuf, mygrid); // every core packs iblank information

    } // mygrid->id_local

    if (id_local > 0) { // if message passing is necessary
      if (mygrid->id_local == id_local) { // non-head ranks send

        // grid point data
        MPI_Send(&count, 1, MPI_INT, 0, id_local, mpi::comm_block); // let its head rank know how many boys I have to send
        MPI_Send(dbuf, count, MPI_DOUBLE, 0, id_local, mpi::comm_block);

        // iblank
        MPI_Send(&(mygrid->num_cells), 1, MPI_INT, 0, id_local, mpi::comm_block); // let its head rank know how many boys I have to send
        MPI_Send(ibuf, mygrid->num_cells, MPI_INT, 0, id_local, mpi::comm_block);

      }
      else if (mygrid->id_local == 0) { // a head rank always receives

        // grid point data
        MPI_Recv(&count, 1, MPI_INT, id_local, id_local, mpi::comm_block, mpi::status);
        DEALLOCATE_1DPTR(dbuf);
        dbuf = new double[count]; // a head rank re-allocates to have the same number of buffer as its sender
        MPI_Recv(dbuf, count, MPI_DOUBLE, id_local, id_local, mpi::comm_block, mpi::status);

        // iblank
        int count_tmp;
        MPI_Recv(&count_tmp, 1, MPI_INT, id_local, id_local, mpi::comm_block, mpi::status);
        DEALLOCATE_1DPTR(ibuf);
        ibuf = new int[count_tmp]; // a head rank re-allocates to have the same number of buffer as its sender
        MPI_Recv(ibuf, count_tmp, MPI_INT, id_local, id_local, mpi::comm_block, mpi::status);

      } // mygrid->id_local
    } // id_local

    if (mpi::irank_block == 0) { // a head rank stores

      int is_in_parent = idxs_in_parent[id_local * DIM_MAX + XI];
      int ie_in_parent = idxe_in_parent[id_local * DIM_MAX + XI];
      int js_in_parent = idxs_in_parent[id_local * DIM_MAX + ETA];
      int je_in_parent = idxe_in_parent[id_local * DIM_MAX + ETA];
      int ks_in_parent = idxs_in_parent[id_local * DIM_MAX + ZETA];
      int ke_in_parent = idxe_in_parent[id_local * DIM_MAX + ZETA];

      int num_data_points = count / num_vars_2write;

      // unpack grid point data
      for (int ivar = 0; ivar < num_vars_2write; ivar++) {

        int num_points = 0;
        for (int k = ks_in_parent; k < ke_in_parent + 1; k++)
          for (int j = js_in_parent; j < je_in_parent + 1; j++)
            for (int i = is_in_parent; i < ie_in_parent + 1; i++) {

              int l0 = block[cur_block].idx1D(i, j, k);

              data_block[num_points_cur_block * ivar + l0] = dbuf[num_data_points * ivar + num_points];
              num_points++;

            } // i
      } // ivar

      // unpack iblank information
      int num_points = 0;
      for (int k = ks_in_parent; k < ke_in_parent + 1; k++)
        for (int j = js_in_parent; j < je_in_parent + 1; j++)
          for (int i = is_in_parent; i < ie_in_parent + 1; i++) {

            int l0 = block[cur_block].idx1D(i, j, k);

            iblank_block[l0] = ibuf[num_points];
            num_points++;

          } // i

    } // mpi::irank_block
  } // id_local
  //
  DEALLOCATE_1DPTR(dbuf);
  DEALLOCATE_1DPTR(ibuf);
  mpi::wait_allothers();

  // fourth, write them all at head ranks
  for (int iblock = 0; iblock < num_blocks; iblock++) { // write in the order of block ID
    if (mygrid->id_local == 0 && iblock == cur_block) { // write only if I am a head rank and only if this is my turn

      std::ofstream ofs;
      ofs.open(cstr_to_constchar(filename), std::ofstream::binary | std::ofstream::app);

      int num_data_points = block[iblock].nXi * block[iblock].nEta * block[iblock].nZeta;

      // write grid point data
      for (int ivar = 0; ivar < num_vars_2write; ivar++) {
        for (int k = 0; k < block[iblock].nZeta; k++)
          for (int j = 0; j < block[iblock].nEta; j++)
            for (int i = 0; i < block[iblock].nXi; i++) {

              int l0 = block[iblock].idx1D(i, j, k);

              ofs.write(reinterpret_cast<const char*>(&data_block[l0 + num_data_points * ivar]), sizeof data_block[0]);

            } // i
      } // ivar

      // write iblank information
      for (int k = 0; k < block[iblock].nZeta; k++)
        for (int j = 0; j < block[iblock].nEta; j++)
          for (int i = 0; i < block[iblock].nXi; i++) {

            int l0 = block[iblock].idx1D(i, j, k);

            ofs.write(reinterpret_cast<const char*>(&iblank_block[l0]), sizeof iblank_block[l0]);

          } // i

      ofs.close();

    } // mygrid->id_local
    mpi::wait_allothers();

  } // iblock

  // now clean up
  if (mpi::irank_block == 0) {// only a head rank in a block allocates

    DEALLOCATE_1DPTR(idxs_in_parent);
    DEALLOCATE_1DPTR(idxe_in_parent);
    DEALLOCATE_1DPTR(data_block);
    DEALLOCATE_1DPTR(iblank_block);

  } // mpi::irank_block
  mpi::wait_allothers();

  return;

} // write_grid_serialIO



void write_solution_serialIO(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, State *mystate, std::string filename_2write) {

  int cur_block, num_blocks, num_points_cur_block;
  int num_grids;
  int num_vars_2write;
  int *idxs_in_parent, *idxe_in_parent;
  double *data_block;
  //
  int count;
  double *dbuf;

  // useful variables
  cur_block = mygrid->id_parent; // a block ID to which this core (or equivalently, grid) belongs
  num_blocks = block[cur_block].num_ourkinds; // number of blocks is shared by every block
  num_grids = mygrid->num_ourkinds;
  num_vars_2write = mystate->num_vars_sol; // mystate->num_vars_sol since solutions are written

  // first, write a header for solution
  write_solution_header(filename_2write, block);

  // second, send block indices to a head rank of each block
  // only a head rank in a block allocates to store
  if (mpi::irank_block == 0) {

    idxs_in_parent = new int[num_grids * DIM_MAX]; // a grid's starting index in the parent block
    idxe_in_parent = new int[num_grids * DIM_MAX]; // a grid's ending index in the parent block

  } // mpi::irank_block
  //
  send_blockindex_2headrank(mygrid, idxs_in_parent, idxe_in_parent);
  mpi::wait_allothers();

  // third, solution information
  // array for local solution data
  count = mygrid->num_cells * num_vars_2write;
  num_points_cur_block = block[cur_block].num_cells;
  dbuf = new double[count];
  if (mpi::irank_block == 0) // only a head rank in a block allocates
    data_block = new double[num_points_cur_block * num_vars_2write];
  //
  // message passing within a block to deliver solution information to a head rank
  for (int id_local = 0; id_local < num_grids; id_local++) {

    if (mygrid->id_local == id_local)
      pack_solution_2send(dbuf, mygrid, mystate); // every core packs solution information

    if (id_local > 0) { // if message passing is necessary
      if (mygrid->id_local == id_local) { // non-head ranks send

        MPI_Send(&count, 1, MPI_INT, 0, id_local, mpi::comm_block); // let its head rank know how many boys I have to send
        MPI_Send(dbuf, count, MPI_DOUBLE, 0, id_local, mpi::comm_block);

      }
      else if (mygrid->id_local == 0) { // a head rank always receives

        MPI_Recv(&count, 1, MPI_INT, id_local, id_local, mpi::comm_block, mpi::status);
        DEALLOCATE_1DPTR(dbuf);
        dbuf = new double[count]; // a head rank re-allocates to have the same number of buffer as its sender
        MPI_Recv(dbuf, count, MPI_DOUBLE, id_local, id_local, mpi::comm_block, mpi::status);

      } // mygrid->id_local
    } // id_local

    if (mpi::irank_block == 0) { // a head rank stores

      int is_in_parent = idxs_in_parent[id_local * DIM_MAX + XI];
      int ie_in_parent = idxe_in_parent[id_local * DIM_MAX + XI];
      int js_in_parent = idxs_in_parent[id_local * DIM_MAX + ETA];
      int je_in_parent = idxe_in_parent[id_local * DIM_MAX + ETA];
      int ks_in_parent = idxs_in_parent[id_local * DIM_MAX + ZETA];
      int ke_in_parent = idxe_in_parent[id_local * DIM_MAX + ZETA];

      int num_data_points = count / num_vars_2write;

      for (int ivar = 0; ivar < num_vars_2write; ivar++) {

        int num_points = 0;
        for (int k = ks_in_parent; k < ke_in_parent + 1; k++)
          for (int j = js_in_parent; j < je_in_parent + 1; j++)
            for (int i = is_in_parent; i < ie_in_parent + 1; i++) {

              int l0 = block[cur_block].idx1D(i, j, k);

              data_block[num_points_cur_block * ivar + l0] = dbuf[num_data_points * ivar + num_points];
              num_points++;

            } // i
      } // ivar
    } // mpi::irank_block
  } // id_local
  //
  DEALLOCATE_1DPTR(dbuf);
  mpi::wait_allothers();

  // fourth, write them all at head ranks
  for (int iblock = 0; iblock < num_blocks; iblock++) { // write in the order of block ID
    if (mygrid->id_local == 0 && iblock == cur_block) { // write only if I am a head rank and only if this is my turn

      std::ofstream ofs;
      ofs.open(cstr_to_constchar(filename_2write), std::ofstream::binary | std::ofstream::app);

      // specific to PLOT3D solution file
      double tau[4];
      tau[0] = 0.0; tau[1] = 1.0; tau[2] = static_cast<double>(mystate->time_step); tau[3] = mystate->time_sol;
      ofs.write(reinterpret_cast<const char*>(&tau), sizeof tau);

      int num_data_points = block[iblock].nXi * block[iblock].nEta * block[iblock].nZeta;
      for (int ivar = 0; ivar < num_vars_2write; ivar++) {
        for (int k = 0; k < block[iblock].nZeta; k++)
          for (int j = 0; j < block[iblock].nEta; j++)
            for (int i = 0; i < block[iblock].nXi; i++) {

              int l0 = block[iblock].idx1D(i, j, k);

              ofs.write(reinterpret_cast<const char*>(&data_block[l0 + num_data_points * ivar]), sizeof data_block[0]);

            } // i
      } // ivar
      ofs.close();

    } // mygrid->id_local
    mpi::wait_allothers();

  } // iblock

  // now clean up
  if (mpi::irank_block == 0) { // only a head rank in a block allocates

    DEALLOCATE_1DPTR(idxs_in_parent);
    DEALLOCATE_1DPTR(idxe_in_parent);
    DEALLOCATE_1DPTR(data_block);

  } // mpi::irank_block
  mpi::wait_allothers();

  return;

} // write_solution_serialIO



void write_function_serialIO(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block, int num_vars_2write, double **func_2write, std::string filename_2write) {

  int cur_block, num_blocks, num_points_cur_block;
  int num_grids;
  int *idxs_in_parent, *idxe_in_parent;
  double *data_block;
  //
  int count;
  double *dbuf;

  // useful variables
  cur_block = mygrid->id_parent; // a block ID to which this core (or equivalently, grid) belongs
  num_blocks = block[cur_block].num_ourkinds; // number of blocks is shared by every block
  num_grids = mygrid->num_ourkinds;

  // first, write a header for function
  write_function_header(filename_2write, block, num_vars_2write);

  // second, send block indices to a head rank of each block
  // only a head rank in a block allocates to store
  if (mpi::irank_block == 0) {

    idxs_in_parent = new int[num_grids * DIM_MAX]; // a grid's starting index in the parent block
    idxe_in_parent = new int[num_grids * DIM_MAX]; // a grid's ending index in the parent block

  } // mpi::irank_block
  //
  send_blockindex_2headrank(mygrid, idxs_in_parent, idxe_in_parent);
  mpi::wait_allothers();

  // third, function information
  // array for local function data
  count = mygrid->num_cells * num_vars_2write;
  num_points_cur_block = block[cur_block].num_cells;
  dbuf = new double[count];
  if (mpi::irank_block == 0) // only a head rank in a block allocates
    data_block = new double[num_points_cur_block * num_vars_2write];
  //
  // message passing within a block to deliver function information to a head rank
  for (int id_local = 0; id_local < num_grids; id_local++) {

    if (mygrid->id_local == id_local)
      pack_function_2send(dbuf, mygrid, num_vars_2write, func_2write);

    if (id_local > 0) { // if message passing is necessary
      if (mygrid->id_local == id_local) { // non-head ranks send

        MPI_Send(&count, 1, MPI_INT, 0, id_local, mpi::comm_block); // let its head rank know how many boys I have to send
        MPI_Send(dbuf, count, MPI_DOUBLE, 0, id_local, mpi::comm_block);

      }
      else if (mygrid->id_local == 0) { // a head rank always receives

        MPI_Recv(&count, 1, MPI_INT, id_local, id_local, mpi::comm_block, mpi::status);
        DEALLOCATE_1DPTR(dbuf);
        dbuf = new double[count]; // a head rank re-allocates to have the same number of buffer as its sender
        MPI_Recv(dbuf, count, MPI_DOUBLE, id_local, id_local, mpi::comm_block, mpi::status);

      } // mygrid->id_local
    } // id_local

    if (mpi::irank_block == 0) { // a head rank stores

      int is_in_parent = idxs_in_parent[id_local * DIM_MAX + XI];
      int ie_in_parent = idxe_in_parent[id_local * DIM_MAX + XI];
      int js_in_parent = idxs_in_parent[id_local * DIM_MAX + ETA];
      int je_in_parent = idxe_in_parent[id_local * DIM_MAX + ETA];
      int ks_in_parent = idxs_in_parent[id_local * DIM_MAX + ZETA];
      int ke_in_parent = idxe_in_parent[id_local * DIM_MAX + ZETA];

      int num_data_points = count / num_vars_2write;

      for (int ivar = 0; ivar < num_vars_2write; ivar++) {

        int num_points = 0;
        for (int k = ks_in_parent; k < ke_in_parent + 1; k++)
          for (int j = js_in_parent; j < je_in_parent + 1; j++)
            for (int i = is_in_parent; i < ie_in_parent + 1; i++) {

              int l0 = block[cur_block].idx1D(i, j, k);

              data_block[num_points_cur_block * ivar + l0] = dbuf[num_data_points * ivar + num_points];
              num_points++;

            } // i
      } // ivar
    } // mpi::irank_block
  } // id_local
  //
  DEALLOCATE_1DPTR(dbuf);
  mpi::wait_allothers();

  // fourth, write them all at head ranks
  for (int iblock = 0; iblock < num_blocks; iblock++) { // write in the order of block ID
    if (mygrid->id_local == 0 && iblock == cur_block) { // write only if I am a head rank and only if this is my turn

      std::ofstream ofs;
      ofs.open(cstr_to_constchar(filename_2write), std::ofstream::binary | std::ofstream::app);

      int num_data_points = block[iblock].nXi * block[iblock].nEta * block[iblock].nZeta;
      for (int ivar = 0; ivar < num_vars_2write; ivar++) {
        for (int k = 0; k < block[iblock].nZeta; k++)
          for (int j = 0; j < block[iblock].nEta; j++)
            for (int i = 0; i < block[iblock].nXi; i++) {

              int l0 = block[iblock].idx1D(i, j, k);
              ofs.write(reinterpret_cast<const char*>(&data_block[l0 + num_data_points * ivar]), sizeof data_block[0]);

            } // i
      } // ivar
      ofs.close();

    } // mygrid->id_local
    mpi::wait_allothers();

  } // iblock

  // now clean up
  if (mpi::irank_block == 0) { // only a head rank in a block allocates

    DEALLOCATE_1DPTR(idxs_in_parent);
    DEALLOCATE_1DPTR(idxe_in_parent);
    DEALLOCATE_1DPTR(data_block);

  } // mpi::irank_block
  mpi::wait_allothers();

  return;

} // write_function_serialIO



void write_solution_namefile(std::string file_varname, State *mystate) {

  write_function_namefile(file_varname, mystate->num_vars_sol, mystate->name_vars);

  return;

} // write_solution_namefile



void write_function_namefile(std::string file_varname, int num_vars, std::string *name_vars) {

  if (mpi::irank == 0) {

    std::ofstream ofs;

    ofs.open(cstr_to_constchar(file_varname), std::ofstream::trunc);
    for (int ivar = 0; ivar < num_vars; ivar++)
      ofs << name_vars[ivar] << std::endl;

    ofs.close();

  } // mpi::irank

  mpi::wait_allothers();

  return;

} // write_function_namefile



std::string *name_metrics() {

  int num_vars = DIM_MAX * DIM_MAX + 2; // 3 by 3 metric tensor, Jacobian, and inverse Jacobian
  std::string *name_vars = new std::string[num_vars];

  name_vars[0] = "XI_X";
  name_vars[1] = "XI_Y";
  name_vars[2] = "XI_Z";
  name_vars[3] = "ETA_X";
  name_vars[4] = "ETA_Y";
  name_vars[5] = "ETA_Z";
  name_vars[6] = "ZETA_X";
  name_vars[7] = "ZETA_Y";
  name_vars[8] = "ZETA_Z";
  name_vars[9] = "J";
  name_vars[10] = "invJ";

  return name_vars;

} // name_metrics

} // plot3d
