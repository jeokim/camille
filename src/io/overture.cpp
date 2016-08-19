#include <iostream>
#include <iomanip>
#include <fstream>

#include "overture.h"

namespace overture {

int num_blocks_in_file = 0;

int num_receiverCells_total = 0;
int *num_receiverCells;

int *iblock_receiver;
int *ijk_receiverCells;
//
int *iblock_donor;
int *ijk_donorCells;
//
int *irank_receiver;
int *irank_donor;

int num_donorBlocks_4myblock = 0;
int *iblock_donor_4myblock;

int num_receiverBlocks_from_myblock = 0;
int *iblock_receiver_from_myblock;

int num_donorRanks_4mygrid = 0;
int *irank_donor_4mygrid;
//
int num_receiverRanks_from_mygrid = 0;
int *irank_receiver_from_mygrid;

int *num_donorCells_4mygrid;
int total_num_donorCells_4mygrid = 0;
//
int *num_receiverCells_from_mygrid;
int total_num_receiverCells_from_mygrid = 0;

int accuracy_interpolation;
int width_stencil[DIM_MAX];
t_InterpStencil *interpStencilExplicit;

// holes
int *num_cells_hole;
int *ijk_holeCells;



void read_overset_NOTSCALING(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  // due to a concern that too many cores attempt to access a single file simultaneously, only head ranks read the file
  // then, the overset information is propagated to grids in each block via message passing

  read_overset_at_headrank_NOTSCALING(myinput, mygrid, block);
  propagate_overset_to_grids_NOTSCALING(myinput, mygrid);

  mpi::wait_allothers();

  return;

} // read_overset_NOTSCALING



void read_overset_at_headrank_NOTSCALING(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  int num_dim;
  std::string filename = myinput->file_overset;
  accuracy_interpolation = myinput->overset_accuracy;

  if (mpi::irank_block == 0) { // only head ranks read

    std::ifstream ifs;

    ifs.open(cstr_to_constchar(filename), std::ifstream::binary | std::ifstream::in);
    if (!ifs.is_open())
      mpi::graceful_exit("The overset-grid file does not exist.");

    // number of dimensions
    ifs.read(reinterpret_cast<char *>(&num_dim), sizeof num_dim);
    assert( num_dim == myinput->num_dim ); // sanity check if grid and overset data have the same dimension
//  std::cout << "number of dimensions: " << std::setw(4) << num_dim << std::endl;

    // number of blocks
    ifs.read(reinterpret_cast<char *>(&num_blocks_in_file), sizeof num_blocks_in_file);
    assert( num_blocks_in_file == block[0].num_ourkinds ); // sanity check if grid and overset data have the same number of blocks
//    std::cout << "num_blocks_in_file: " << std::setw(4) << num_blocks_in_file << std::endl;

    // number of cells per block
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      int num_cells_in_this_block;
      ifs.read(reinterpret_cast<char *>(&num_cells_in_this_block), sizeof num_cells_in_this_block);
      assert( num_cells_in_this_block == block[iblock].num_cells ); // sanity check if grid and overset data have the same cells per block
//      std::cout << "Block ID: " << std::setw(4) << iblock << " , number of cells: " << std::setw(8) << num_cells_in_this_block << std::endl;

    } // iblock

    // number of interpolation points (receiver points)
    num_receiverCells_total = 0;
    num_receiverCells = new int[num_blocks_in_file];
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      ifs.read(reinterpret_cast<char *>(&num_receiverCells[iblock]), sizeof num_receiverCells[iblock]);
//      std::cout << "Block ID: " << std::setw(4) << iblock << ", global rank: " << std::setw(4) << mpi::irank << " , number of receivers: " << std::setw(8) << num_receiverCells[iblock] << std::endl;

      num_receiverCells_total += num_receiverCells[iblock];

    } // iblock

    // receiver-block indices
    iblock_receiver = new int[num_receiverCells_total];
    for (int counter = 0; counter < num_receiverCells_total; counter++)
      ifs.read(reinterpret_cast<char *>(&iblock_receiver[counter]), sizeof iblock_receiver[counter]);

    // block-level indices for the interpolation points (receiver points)
    ijk_receiverCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_receiverCells[counter] = 0;
    //
    int counter1 = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
      for (int ircv = 0; ircv < num_receiverCells[iblock]; ircv++) {
        for (int idir = XI; idir < DIM_MAX; idir++) {

          if (idir < num_dim) // bypass the file-reading in a direction higher than the number of dimension of interest
            ifs.read(reinterpret_cast<char *>(&ijk_receiverCells[counter1]), sizeof ijk_receiverCells[counter1]);
          counter1++;

        } // idir
      } // ircv
    } // iblock
    //
//    counter1 = 0;
//    for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
//      for (int ircv = 0; ircv < num_receiverCells[iblock]; ircv++)
//        for (int idir = XI; idir < DIM_MAX; idir++) {
//          if (mpi::irank == 0)
//            std::cout << "Block: "  << std::setw(4) << iblock
//                      << ", ircv: " << std::setw(8) << ircv
//                      << ", idir: " << std::setw(4) << idir
//                      << ", ijk: "  << std::setw(8) << ijk_receiverCells[counter1] << std::endl;
//          counter1++;
//        } // idir
    assert( counter1 == num_receiverCells_total * DIM_MAX);

    // donor-block indices
    iblock_donor = new int[num_receiverCells_total];
    for (int counter = 0; counter < num_receiverCells_total; counter++) {

      ifs.read(reinterpret_cast<char *>(&iblock_donor[counter]), sizeof iblock_donor[counter]);
//      if (mpi::irank == 0)
//        std::cout << "Receiver-cell index: " << counter << ", iblock_donor: " << iblock_donor[counter] << std::endl;

    } // counter
    //
    identify_my_donor_blocks(myinput, mygrid);

    // block-level indices for the donor cells (donor points)
    ijk_donorCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_donorCells[counter] = 0;
    int counter2 = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
      for (int idonor = 0; idonor < num_receiverCells[iblock]; idonor++) { // there is a one-to-one correspondence between receivers and donors
        for (int idir = XI; idir < DIM_MAX; idir++) {

          if (idir < num_dim) // bypass the file-reading in a direction higher than the number of dimension of interest
            ifs.read(reinterpret_cast<char *>(&ijk_donorCells[counter2]), sizeof ijk_donorCells[counter2]);
          counter2++;

        } // idir
      } // idonor
    } // iblock
    //
//    counter2 = 0;
//    for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
//      for (int idonor = 0; idonor < num_receiverCells[iblock]; idonor++)
//        for (int idir = XI; idir < DIM_MAX; idir++) {
//          if (mpi::irank == 0)
//            std::cout << "Block: "  << std::setw(4) << iblock
//                      << ", idonor: " << std::setw(8) << idonor
//                      << ", idir: " << std::setw(4) << idir
//                      << ", ijk: "  << std::setw(8) << ijk_donorCells[counter2] << std::endl;
//          counter2++;
//        } // idir
    assert( counter2 == num_receiverCells_total * DIM_MAX);

    identify_my_receiver_blocks(myinput, mygrid);

    // interpolation weights
    interpStencilExplicit = new t_InterpStencil[num_receiverCells_total]; // for now, store the stencil coefficients for every donor point
    //
    int icell_base = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      // size of interpolation stencils
      for (int idir = XI; idir < DIM_MAX; idir++)
        ifs.read(reinterpret_cast<char *>(&width_stencil[idir]), sizeof width_stencil[idir]);
      for (int idir = XI; idir < myinput->num_dim; idir++)
        assert( width_stencil[idir] - 1 == accuracy_interpolation ); // consistency check between the user input (myinput->overset_accuracy) and what's in the Overture file
      //
//      std::cout << "Block: " << std::setw(4) << iblock << ", width_XI"   << std::setw(4) << width_stencil[XI]
//                                                       << ", width_ETA"  << std::setw(4) << width_stencil[ETA]
//                                                       << ", width_ZETA" << std::setw(4) << width_stencil[ZETA] << std::endl;

      int stencilSize = math_algebra::product_array_integer(width_stencil, DIM_MAX);

      // interpolation coefficients
		  for (int icell = icell_base; icell < icell_base + num_receiverCells[iblock]; icell++) {

        interpStencilExplicit[icell].coeff = new double[stencilSize];

        int counter_stencil = 0;
	      for (int w3 = 0; w3 < width_stencil[ZETA]; w3++)
		      for (int w2 = 0; w2 < width_stencil[ETA]; w2++)
		        for (int w1 = 0; w1 < width_stencil[XI]; w1++) {

              double coeff = 0.0;
              ifs.read(reinterpret_cast<char *>(&coeff), sizeof coeff);
              interpStencilExplicit[icell].coeff[counter_stencil] = coeff;
//              if (mpi::irank == 0)
//                std::cout << "Block: "     << std::setw(4) << iblock
//                          << ", icell: "   << std::setw(4) << icell - icell_base
//                          << ", counter: " << std::setw(4) << counter_stencil
//                          << ", coeff: "   << interpStencilExplicit[icell].coeff[counter_stencil] << std::endl;
              counter_stencil++;

            } // w1
      } // icell
      icell_base += num_receiverCells[iblock];

    } // iblock

    // how many hole cells we have per block
    num_cells_hole = new int[num_blocks_in_file];
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      ifs.read(reinterpret_cast<char *>(&num_cells_hole[iblock]), sizeof num_cells_hole[iblock]);
//      std::cout << "Rank: " << std::setw(4) << mpi::irank
//                << ", block: " << std::setw(4) << iblock
//                << ", num_cells_hole: " << std::setw(8) << num_cells_hole[iblock] << std::endl;

    } // iblock
    int num_cells_hole_total = math_algebra::sum_array_integer(num_cells_hole, num_blocks_in_file);

    // block-level indices for masked cells
    ijk_holeCells = new int[num_cells_hole_total * DIM_MAX];
    int num_cells;
    int counter = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
      for (int ihole = 0; ihole < num_cells_hole[iblock]; ihole++) {
        for (int idir = XI; idir < DIM_MAX; idir++) {

          ifs.read(reinterpret_cast<char *>(&num_cells), sizeof num_cells);
          ijk_holeCells[idir + counter * DIM_MAX] = num_cells;

        } // idir
        counter++;

      } // ihole
    } // iblock

    ifs.close();

  } // mpi::irank_block

  mpi::wait_allothers();

  return;

} // read_overset_at_headrank_NOTSCALING



void propagate_overset_to_grids_NOTSCALING(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // note that every communication occurs within each block

  MPI_Bcast(&num_blocks_in_file, 1, MPI_INT, 0, mpi::comm_block);

  MPI_Bcast(&num_receiverCells_total, 1, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    num_receiverCells = new int[num_blocks_in_file];
  MPI_Bcast(num_receiverCells, num_blocks_in_file, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    iblock_receiver = new int[num_receiverCells_total];
  MPI_Bcast(iblock_receiver, num_receiverCells_total, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    ijk_receiverCells = new int[num_receiverCells_total * DIM_MAX];
  MPI_Bcast(ijk_receiverCells, num_receiverCells_total * DIM_MAX, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    iblock_donor = new int[num_receiverCells_total];
  MPI_Bcast(iblock_donor, num_receiverCells_total, MPI_INT, 0, mpi::comm_block);

  MPI_Bcast(&num_donorBlocks_4myblock, 1, MPI_INT, 0, mpi::comm_block);
  if (mpi::irank_block > 0)
    iblock_donor_4myblock = new int[num_donorBlocks_4myblock];
  MPI_Bcast(iblock_donor_4myblock, num_donorBlocks_4myblock, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    ijk_donorCells = new int[num_receiverCells_total * DIM_MAX];
  MPI_Bcast(ijk_donorCells, num_receiverCells_total * DIM_MAX, MPI_INT, 0, mpi::comm_block);

  MPI_Bcast(&num_receiverBlocks_from_myblock, 1, MPI_INT, 0, mpi::comm_block);
  if (mpi::irank_block > 0)
    iblock_receiver_from_myblock = new int[num_receiverBlocks_from_myblock];
  MPI_Bcast(iblock_receiver_from_myblock, num_receiverBlocks_from_myblock, MPI_INT, 0, mpi::comm_block);

  MPI_Bcast(width_stencil, DIM_MAX, MPI_INT, 0, mpi::comm_block);
  int stencilSize = math_algebra::product_array_integer(width_stencil, DIM_MAX);

  if (mpi::irank_block > 0) {

    interpStencilExplicit = new t_InterpStencil[num_receiverCells_total];

  	for (int icell = 0; icell < num_receiverCells_total; icell++)
      interpStencilExplicit[icell].coeff = new double[stencilSize];

  } // mpi::irank_block > 0
  for (int icell = 0; icell < num_receiverCells_total; icell++)
    MPI_Bcast(interpStencilExplicit[icell].coeff, stencilSize, MPI_DOUBLE, 0, mpi::comm_block);
  //
//  if (mpi::irank == 0)
//    for (int icell = 0; icell < num_receiverCells_total; icell++)
//      for (int istencil = 0; istencil < stencilSize; istencil++)
//        std::cout << "Rank: " << std::setw(4) << mpi::irank << ", icell: " << std::setw(4) << icell << ", istencil: " << std::setw(4) << istencil
//                  << ", coeff: " << std::scientific << interpStencilExplicit[icell].coeff[istencil] << std::endl;

  if (mpi::irank_block > 0)
    num_cells_hole = new int[num_blocks_in_file];
  MPI_Bcast(num_cells_hole, num_blocks_in_file, MPI_INT, 0, mpi::comm_block);

  int num_cells_hole_total = math_algebra::sum_array_integer(num_cells_hole, num_blocks_in_file);
  if (mpi::irank_block > 0)
    ijk_holeCells = new int[num_cells_hole_total * DIM_MAX];
  MPI_Bcast(ijk_holeCells, num_cells_hole_total * DIM_MAX, MPI_INT, 0, mpi::comm_block);

  return;

} // propagate_overset_to_grids_NOTSCALING



void identify_my_donor_blocks(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  if (mpi::irank_block != 0)
    return;

  int cur_block = mygrid->id_parent;
  int counter_base;
  int *if_my_donor_blocks = new int[num_blocks_in_file];

  for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
    if_my_donor_blocks[iblock] = FALSE;

  counter_base = 0;
  for (int iblock = 0; iblock < cur_block; iblock++)
    counter_base += num_receiverCells[iblock];

  for (int counter = 0; counter < num_receiverCells[cur_block]; counter++) { // check my portion of iblock_donor

    int index_block = iblock_donor[counter_base + counter];
    if_my_donor_blocks[index_block] = TRUE;

  } // counter

  // count how many blocks are sending data to me
  for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
    if (if_my_donor_blocks[iblock] == TRUE)
      num_donorBlocks_4myblock++;

  // store the indices of the blocks that I receive for later message passing
  iblock_donor_4myblock = new int[num_donorBlocks_4myblock];
  int num_donor_blocks = 0;
  for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
    if (if_my_donor_blocks[iblock] == TRUE) {

      iblock_donor_4myblock[num_donor_blocks] = iblock;
      num_donor_blocks++;

    } // if_my_donor_blocks[iblock]
//  std::cout << "Rank: " << mpi::irank << ", num_donor_blocks: " << num_donor_blocks << std::endl;
//  for (int iblock = 0; iblock < num_donorBlocks_4myblock; iblock++)
//    std::cout << "Rank: " << mpi::irank << ", iblock: " << iblock << ", iblock_donor_4myblock: " << iblock_donor_4myblock[iblock] << std::endl;
  assert( num_donor_blocks == num_donorBlocks_4myblock );

  DEALLOCATE_1DPTR(if_my_donor_blocks);

  return;

} // identify_my_donor_blocks



void identify_my_receiver_blocks(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  if (mpi::irank_block != 0)
    return;

  int cur_block = mygrid->id_parent;
  int *if_my_receiver_blocks = new int[num_blocks_in_file];

  for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
    if_my_receiver_blocks[iblock] = FALSE;

  int counter = 0;
  for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
    for (int ircv = 0; ircv < num_receiverCells[iblock]; ircv++) {

      if (iblock_donor[counter] == cur_block) // if I am a donor block, the other block must be a receiver
        if_my_receiver_blocks[iblock] = TRUE;

      counter++;

    } // ircv
  } // iblock

  // count how many blocks I have to send data to
  for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
    if (if_my_receiver_blocks[iblock] == TRUE)
      num_receiverBlocks_from_myblock++;

  // store the indices of the blocks that I donate for later message passing
  iblock_receiver_from_myblock = new int[num_receiverBlocks_from_myblock];
  int num_receiving_blocks = 0;
  for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
    if (if_my_receiver_blocks[iblock] == TRUE) {

      iblock_receiver_from_myblock[num_receiving_blocks] = iblock;
      num_receiving_blocks++;

    } // if_my_receiver_blocks[iblock]
//  std::cout << "Rank: " << mpi::irank << ", num_receiving_blocks: " << num_receiving_blocks << std::endl;
//  for (int iblock = 0; iblock < num_receiverBlocks_from_myblock; iblock++)
//    std::cout << "Rank: " << mpi::irank << ", iblock: " << iblock << ", iblock_receiver_from_myblock: " << iblock_receiver_from_myblock[iblock] << std::endl;
  assert( num_receiving_blocks == num_receiverBlocks_from_myblock );

  DEALLOCATE_1DPTR(if_my_receiver_blocks);

  return;

} // identify_my_receiver_blocks



void read_overset(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  // due to a concern that too many cores attempt to access a single file simultaneously, only head ranks read the file
  // then, the overset information is propagated to grids in each block via message passing

  read_overset_at_headrank(myinput, mygrid, block);
  propagate_overset_to_grids(myinput, mygrid);

  convert_blockLevelInfo2GridLevel(myinput, mygrid);
  identify_my_donor_ranks(myinput, mygrid);
  identify_my_receiver_ranks(myinput, mygrid);

  mpi::wait_allothers();

  return;

} // read_overset



void read_overset_at_headrank(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  int num_dim;
  std::string filename = myinput->file_overset;
  accuracy_interpolation = myinput->overset_accuracy;

  if (mpi::irank_block == 0) { // only head ranks read

    std::ifstream ifs;

    ifs.open(cstr_to_constchar(filename), std::ifstream::binary | std::ifstream::in);
    if (!ifs.is_open())
      mpi::graceful_exit("The overset-grid file does not exist.");

    // number of dimensions
    ifs.read(reinterpret_cast<char *>(&num_dim), sizeof num_dim);
    assert( num_dim == myinput->num_dim ); // sanity check if grid and overset data have the same dimension
//  std::cout << "number of dimensions: " << std::setw(4) << num_dim << std::endl;

    // number of blocks
    ifs.read(reinterpret_cast<char *>(&num_blocks_in_file), sizeof num_blocks_in_file);
    assert( num_blocks_in_file == block[0].num_ourkinds ); // sanity check if grid and overset data have the same number of blocks
//    std::cout << "num_blocks_in_file: " << std::setw(4) << num_blocks_in_file << std::endl;

    // number of cells per block
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      int num_cells_in_this_block;
      ifs.read(reinterpret_cast<char *>(&num_cells_in_this_block), sizeof num_cells_in_this_block);
      assert( num_cells_in_this_block == block[iblock].num_cells ); // sanity check if grid and overset data have the same cells per block
//      std::cout << "Block ID: " << std::setw(4) << iblock << " , number of cells: " << std::setw(8) << num_cells_in_this_block << std::endl;

    } // iblock

    // number of interpolation points (receiver points)
    num_receiverCells_total = 0;
    num_receiverCells = new int[num_blocks_in_file];
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      ifs.read(reinterpret_cast<char *>(&num_receiverCells[iblock]), sizeof num_receiverCells[iblock]);
//      std::cout << "Block ID: " << std::setw(4) << iblock << ", global rank: " << std::setw(4) << mpi::irank << " , number of receivers: " << std::setw(8) << num_receiverCells[iblock] << std::endl;

      num_receiverCells_total += num_receiverCells[iblock];

    } // iblock

    // receiver-block indices
    iblock_receiver = new int[num_receiverCells_total];
    for (int counter = 0; counter < num_receiverCells_total; counter++)
      ifs.read(reinterpret_cast<char *>(&iblock_receiver[counter]), sizeof iblock_receiver[counter]);

    // block-level indices for the interpolation points (receiver points)
    ijk_receiverCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_receiverCells[counter] = 0;
    //
    int counter1 = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
      for (int ircv = 0; ircv < num_receiverCells[iblock]; ircv++) {
        for (int idir = XI; idir < DIM_MAX; idir++) {

          if (idir < num_dim) // bypass the file-reading in a direction higher than the number of dimension of interest
            ifs.read(reinterpret_cast<char *>(&ijk_receiverCells[counter1]), sizeof ijk_receiverCells[counter1]);
          counter1++;

        } // idir
      } // ircv
    } // iblock
    //
//    counter1 = 0;
//    for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
//      for (int ircv = 0; ircv < num_receiverCells[iblock]; ircv++)
//        for (int idir = XI; idir < DIM_MAX; idir++) {
//          if (mpi::irank == 0)
//            std::cout << "Block: "  << std::setw(4) << iblock
//                      << ", ircv: " << std::setw(8) << ircv
//                      << ", idir: " << std::setw(4) << idir
//                      << ", ijk: "  << std::setw(8) << ijk_receiverCells[counter1] << std::endl;
//          counter1++;
//        } // idir
    assert( counter1 == num_receiverCells_total * DIM_MAX);

    // donor-block indices
    iblock_donor = new int[num_receiverCells_total];
    for (int counter = 0; counter < num_receiverCells_total; counter++) {

      ifs.read(reinterpret_cast<char *>(&iblock_donor[counter]), sizeof iblock_donor[counter]);
//      if (mpi::irank == 0)
//        std::cout << "Receiver-cell index: " << counter << ", iblock_donor: " << iblock_donor[counter] << std::endl;

    } // counter
    //
//    identify_my_donor_blocks(myinput, mygrid);

    // block-level indices for the donor cells (donor points)
    ijk_donorCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_donorCells[counter] = 0;
    int counter2 = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
      for (int idonor = 0; idonor < num_receiverCells[iblock]; idonor++) { // there is a one-to-one correspondence between receivers and donors
        for (int idir = XI; idir < DIM_MAX; idir++) {

          if (idir < num_dim) // bypass the file-reading in a direction higher than the number of dimension of interest
            ifs.read(reinterpret_cast<char *>(&ijk_donorCells[counter2]), sizeof ijk_donorCells[counter2]);
          counter2++;

        } // idir
      } // idonor
    } // iblock
    //
//    counter2 = 0;
//    for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
//      for (int idonor = 0; idonor < num_receiverCells[iblock]; idonor++)
//        for (int idir = XI; idir < DIM_MAX; idir++) {
//          if (mpi::irank == 0)
//            std::cout << "Block: "  << std::setw(4) << iblock
//                      << ", idonor: " << std::setw(8) << idonor
//                      << ", idir: " << std::setw(4) << idir
//                      << ", ijk: "  << std::setw(8) << ijk_donorCells[counter2] << std::endl;
//          counter2++;
//        } // idir
    assert( counter2 == num_receiverCells_total * DIM_MAX);

//    identify_my_receiver_blocks(myinput, mygrid);

    // interpolation weights
    interpStencilExplicit = new t_InterpStencil[num_receiverCells_total]; // for now, store the stencil coefficients for every donor point
    //
    int icell_base = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      // size of interpolation stencils
      for (int idir = XI; idir < DIM_MAX; idir++)
        ifs.read(reinterpret_cast<char *>(&width_stencil[idir]), sizeof width_stencil[idir]);
      for (int idir = XI; idir < myinput->num_dim; idir++)
        assert( width_stencil[idir] - 1 == accuracy_interpolation ); // consistency check between the user input (myinput->overset_accuracy) and what's in the Overture file
      //
//      std::cout << "Block: " << std::setw(4) << iblock << ", width_XI"   << std::setw(4) << width_stencil[XI]
//                                                       << ", width_ETA"  << std::setw(4) << width_stencil[ETA]
//                                                       << ", width_ZETA" << std::setw(4) << width_stencil[ZETA] << std::endl;

      int stencilSize = math_algebra::product_array_integer(width_stencil, DIM_MAX);

      // interpolation coefficients
		  for (int icell = icell_base; icell < icell_base + num_receiverCells[iblock]; icell++) {

        interpStencilExplicit[icell].coeff = new double[stencilSize];

        int counter_stencil = 0;
	      for (int w3 = 0; w3 < width_stencil[ZETA]; w3++)
		      for (int w2 = 0; w2 < width_stencil[ETA]; w2++)
		        for (int w1 = 0; w1 < width_stencil[XI]; w1++) {

              double coeff = 0.0;
              ifs.read(reinterpret_cast<char *>(&coeff), sizeof coeff);
              interpStencilExplicit[icell].coeff[counter_stencil] = coeff;
//              if (mpi::irank == 0)
//                std::cout << "Block: "     << std::setw(4) << iblock
//                          << ", icell: "   << std::setw(4) << icell - icell_base
//                          << ", counter: " << std::setw(4) << counter_stencil
//                          << ", coeff: "   << interpStencilExplicit[icell].coeff[counter_stencil] << std::endl;
              counter_stencil++;

            } // w1
      } // icell
      icell_base += num_receiverCells[iblock];

    } // iblock

    // how many hole cells we have per block
    num_cells_hole = new int[num_blocks_in_file];
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {

      ifs.read(reinterpret_cast<char *>(&num_cells_hole[iblock]), sizeof num_cells_hole[iblock]);
//      std::cout << "Rank: " << std::setw(4) << mpi::irank
//                << ", block: " << std::setw(4) << iblock
//                << ", num_cells_hole: " << std::setw(8) << num_cells_hole[iblock] << std::endl;

    } // iblock
    int num_cells_hole_total = math_algebra::sum_array_integer(num_cells_hole, num_blocks_in_file);

    // block-level indices for masked cells
    ijk_holeCells = new int[num_cells_hole_total * DIM_MAX];
    int num_cells;
    int counter = 0;
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
      for (int ihole = 0; ihole < num_cells_hole[iblock]; ihole++) {
        for (int idir = XI; idir < DIM_MAX; idir++) {

          ifs.read(reinterpret_cast<char *>(&num_cells), sizeof num_cells);
          ijk_holeCells[idir + counter * DIM_MAX] = num_cells;

        } // idir
        counter++;

      } // ihole
    } // iblock

    ifs.close();

  } // mpi::irank_block

  mpi::wait_allothers();

  return;

} // read_overset_at_headrank



void propagate_overset_to_grids(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // note that every communication occurs within each block

  MPI_Bcast(&num_blocks_in_file, 1, MPI_INT, 0, mpi::comm_block);

  MPI_Bcast(&num_receiverCells_total, 1, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    num_receiverCells = new int[num_blocks_in_file];
  MPI_Bcast(num_receiverCells, num_blocks_in_file, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    iblock_receiver = new int[num_receiverCells_total];
  MPI_Bcast(iblock_receiver, num_receiverCells_total, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    ijk_receiverCells = new int[num_receiverCells_total * DIM_MAX];
  MPI_Bcast(ijk_receiverCells, num_receiverCells_total * DIM_MAX, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    iblock_donor = new int[num_receiverCells_total];
  MPI_Bcast(iblock_donor, num_receiverCells_total, MPI_INT, 0, mpi::comm_block);

//  MPI_Bcast(&num_donorBlocks_4myblock, 1, MPI_INT, 0, mpi::comm_block);
//  if (mpi::irank_block > 0)
//    iblock_donor_4myblock = new int[num_donorBlocks_4myblock];
//  MPI_Bcast(iblock_donor_4myblock, num_donorBlocks_4myblock, MPI_INT, 0, mpi::comm_block);

  if (mpi::irank_block > 0)
    ijk_donorCells = new int[num_receiverCells_total * DIM_MAX];
  MPI_Bcast(ijk_donorCells, num_receiverCells_total * DIM_MAX, MPI_INT, 0, mpi::comm_block);

//  MPI_Bcast(&num_receiverBlocks_from_myblock, 1, MPI_INT, 0, mpi::comm_block);
//  if (mpi::irank_block > 0)
//    iblock_receiver_from_myblock = new int[num_receiverBlocks_from_myblock];
//  MPI_Bcast(iblock_receiver_from_myblock, num_receiverBlocks_from_myblock, MPI_INT, 0, mpi::comm_block);

  MPI_Bcast(width_stencil, DIM_MAX, MPI_INT, 0, mpi::comm_block);
  int stencilSize = math_algebra::product_array_integer(width_stencil, DIM_MAX);

  if (mpi::irank_block > 0) {

    interpStencilExplicit = new t_InterpStencil[num_receiverCells_total];

  	for (int icell = 0; icell < num_receiverCells_total; icell++)
      interpStencilExplicit[icell].coeff = new double[stencilSize];

  } // mpi::irank_block > 0
  for (int icell = 0; icell < num_receiverCells_total; icell++)
    MPI_Bcast(interpStencilExplicit[icell].coeff, stencilSize, MPI_DOUBLE, 0, mpi::comm_block);
  //
//  if (mpi::irank == 0)
//    for (int icell = 0; icell < num_receiverCells_total; icell++)
//      for (int istencil = 0; istencil < stencilSize; istencil++)
//        std::cout << "Rank: " << std::setw(4) << mpi::irank << ", icell: " << std::setw(4) << icell << ", istencil: " << std::setw(4) << istencil
//                  << ", coeff: " << std::scientific << interpStencilExplicit[icell].coeff[istencil] << std::endl;

  if (mpi::irank_block > 0)
    num_cells_hole = new int[num_blocks_in_file];
  MPI_Bcast(num_cells_hole, num_blocks_in_file, MPI_INT, 0, mpi::comm_block);

  int num_cells_hole_total = math_algebra::sum_array_integer(num_cells_hole, num_blocks_in_file);
  if (mpi::irank_block > 0)
    ijk_holeCells = new int[num_cells_hole_total * DIM_MAX];
  MPI_Bcast(ijk_holeCells, num_cells_hole_total * DIM_MAX, MPI_INT, 0, mpi::comm_block);

  return;

} // propagate_overset_to_grids



void convert_blockLevelInfo2GridLevel(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // At this point, every core (even if it does not have interpolation points nor 
  // donor points) knows all the overset-interpolation information.  In the initial 
  // implementation (see those with OVERTURE_NOTSCALING), cores within each block 
  // interpolated, the interpolated solutions are merged into the head rank, then 
  // only head ranks did message passing.  This did not scale as the number of cores 
  // and the number of interpolation points increase.  For scalable interpolation on-the-fly, 
  // some of the information should be converted from block level to grid level, especially, 
  //   iblock_receiver: block IDs for interpolation points
  //     -> irank_receiver: core ranks for interpolation points
  //   iblock_donor: block IDs for donor stencils
  //     -> irank_donor: core ranks for donor stencils

  int *buf_sum = new int[num_receiverCells_total];
  int howsitgoing;

  // get irank_receiver out of iblock_receiver and ijk_receiverCells
  irank_receiver = new int[num_receiverCells_total];
  for (int counter = 0; counter < num_receiverCells_total; counter++) {

    irank_receiver[counter] = 0;

    // check if this receiver point is in my block
    if (iblock_receiver[counter] != mygrid->id_parent)
      continue;

    // check if this receiver point is in my grid
    int found = TRUE;
    int idir = XI;
    while ( found == TRUE && idir < DIM_MAX) {

      int counter2 = DIM_MAX * counter + idir;
      if (ijk_receiverCells[counter2] >= mygrid->is_in_parent[idir] && 
          ijk_receiverCells[counter2] <= mygrid->ie_in_parent[idir]) {

        found = TRUE;
        idir++;

      } // ijk_receiverCells[counter2]
      else
        found = FALSE;

    } // found

    if (found == TRUE)
      irank_receiver[counter] = mpi::irank + 1; // temporarily increased by 1

  } // counter

  // sum irank_receiver over all the cores so that the other cores know my contribution
  for (int counter = 0; counter < num_receiverCells_total; counter++)
    buf_sum[counter] = 0;
  //
  MPI_Allreduce(irank_receiver, buf_sum, num_receiverCells_total, MPI_INT, MPI_SUM, mpi::comm_region);
  //
  howsitgoing = OK;
  for (int counter = 0; counter < num_receiverCells_total; counter++) {

    irank_receiver[counter] = buf_sum[counter] - 1; // now decreased by 1
    if (irank_receiver[counter] < 0 || irank_receiver[counter] >= mpi::nprocs)
      howsitgoing = NOT_OK;

  } // counter
  if (howsitgoing != OK)
    mpi::graceful_exit("There is at least one receiver point which does not belong to any core.");
  //
  mpi::wait_allothers();

  // get irank_donor out of iblock_donor and ijk_donorCells
  irank_donor = new int[num_receiverCells_total];
  for (int counter = 0; counter < num_receiverCells_total; counter++) {

    irank_donor[counter] = 0;

    // check if this donor point is in my block
    if (iblock_donor[counter] != mygrid->id_parent)
      continue;

    // check if this donor point is in my grid
    int found = TRUE;
    int idir = XI;
    while ( found == TRUE && idir < DIM_MAX) {

      int counter2 = DIM_MAX * counter + idir;
      if (ijk_donorCells[counter2] >= mygrid->is_in_parent[idir] && 
          ijk_donorCells[counter2] <= mygrid->ie_in_parent[idir]) {

        found = TRUE;
        idir++;

      } // ijk_donorCells[counter2]
      else
        found = FALSE;

    } // found

    if (found == TRUE)
      irank_donor[counter] = mpi::irank + 1; // temporarily increased by 1

  } // counter

  // sum irank_donor over all the cores so that the other cores know my contribution
  for (int counter = 0; counter < num_receiverCells_total; counter++)
    buf_sum[counter] = 0;
  //
  MPI_Allreduce(irank_donor, buf_sum, num_receiverCells_total, MPI_INT, MPI_SUM, mpi::comm_region);
  //
  howsitgoing = OK;
  for (int counter = 0; counter < num_receiverCells_total; counter++) {

    irank_donor[counter] = buf_sum[counter] - 1; // now decreased by 1
    if (irank_donor[counter] < 0 || irank_donor[counter] >= mpi::nprocs)
      howsitgoing = NOT_OK;

  } // counter
  if (howsitgoing != OK)
    mpi::graceful_exit("There is at least one donor stencil which does not belong to any core.");
  //
  mpi::wait_allothers();

  DEALLOCATE_1DPTR(buf_sum);

  return;

} // convert_blockLevelInfo2GridLevel



void identify_my_donor_ranks(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  int *if_my_donor_ranks = new int[mpi::nprocs];

  for (int irank = 0; irank < mpi::nprocs; irank++)
    if_my_donor_ranks[irank] = FALSE; // by default, no rank is my donor

  for (int counter = 0; counter < num_receiverCells_total; counter++)
    if (irank_receiver[counter] == mpi::irank) // if this receiver point belongs to me
      if_my_donor_ranks[irank_donor[counter]] = TRUE;

  // count how many cores are sending data to me
  for (int irank = 0; irank < mpi::nprocs; irank++)
    if (if_my_donor_ranks[irank] == TRUE)
      num_donorRanks_4mygrid++;

  // store the ranks from which I receive for later message passing
  irank_donor_4mygrid = new int[num_donorRanks_4mygrid];
  int num_donor_ranks = 0;
  for (int irank = 0; irank < mpi::nprocs; irank++)
    if (if_my_donor_ranks[irank] == TRUE) {

      irank_donor_4mygrid[num_donor_ranks] = irank;
      num_donor_ranks++;

    } // if_my_donor_ranks[irank]
  assert( num_donor_ranks == num_donorRanks_4mygrid );

  // store the number of donor cells in my donor ranks
  num_donorCells_4mygrid = new int[num_donorRanks_4mygrid];
  for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
    num_donorCells_4mygrid[irank] = 0;
  //
  for (int counter = 0; counter < num_receiverCells_total; counter++)
    if (irank_receiver[counter] == mpi::irank) { // if this receiver point belongs to me

      // the below routine is equivalent to finding num_donor_ranks given irank for 
      // irank_donor_4mygrid[num_donor_ranks] = irank
      int this_index = NONE;
      for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
        if (irank_donor_4mygrid[irank] == irank_donor[counter])
          this_index = irank;
      assert( this_index != NONE);

      num_donorCells_4mygrid[this_index]++;

    } // irank_receiver[counter]
  //
  for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
    assert( num_donorCells_4mygrid[irank] > 0 );

  for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
    total_num_donorCells_4mygrid += num_donorCells_4mygrid[irank];

  DEALLOCATE_1DPTR(if_my_donor_ranks);

  return;

} // identify_my_donor_ranks



void identify_my_receiver_ranks(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  int *if_my_receiver_ranks = new int[mpi::nprocs];

  for (int irank = 0; irank < mpi::nprocs; irank++)
    if_my_receiver_ranks[irank] = FALSE; // by default, no rank is my receiver

  for (int counter = 0; counter < num_receiverCells_total; counter++)
    if (irank_donor[counter] == mpi::irank) // if this donor point belongs to me
      if_my_receiver_ranks[irank_receiver[counter]] = TRUE;

  // count how many cores I have to send data to
  for (int irank = 0; irank < mpi::nprocs; irank++)
    if (if_my_receiver_ranks[irank] == TRUE)
      num_receiverRanks_from_mygrid++;

  // store the ranks to which I donate for later message passing
  irank_receiver_from_mygrid = new int[num_receiverRanks_from_mygrid];
  int num_receiver_ranks = 0;
  for (int irank = 0; irank < mpi::nprocs; irank++)
    if (if_my_receiver_ranks[irank] == TRUE) {

      irank_receiver_from_mygrid[num_receiver_ranks] = irank;
      num_receiver_ranks++;

    } // if_my_receiver_ranks[irank]
  assert( num_receiver_ranks == num_receiverRanks_from_mygrid );

  // store the number of receiver cells in my receiver ranks
  num_receiverCells_from_mygrid = new int[num_receiverRanks_from_mygrid];
  for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
    num_receiverCells_from_mygrid[irank] = 0;
  //
  for (int counter = 0; counter < num_receiverCells_total; counter++)
    if (irank_donor[counter] == mpi::irank) { // if this donor point point belongs to me

      // the below routine is equivalent to finding num_receiver_ranks given irank for 
      // irank_receiver_from_mygrid[num_receiver_ranks] = irank
      int this_index = NONE;
      for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
        if (irank_receiver_from_mygrid[irank] == irank_receiver[counter])
          this_index = irank;
      assert( this_index != NONE);

      num_receiverCells_from_mygrid[this_index]++;

    } // irank_donor[counter]
  //
  for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
    assert( num_receiverCells_from_mygrid[irank] > 0 );

  for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
    total_num_receiverCells_from_mygrid += num_receiverCells_from_mygrid[irank];

  DEALLOCATE_1DPTR(if_my_receiver_ranks);

  return;

} // identify_my_receiver_ranks



void cleanup_overset_temporarydata(UserInput *myinput) {

  DEALLOCATE_1DPTR(interpStencilExplicit);

  num_blocks_in_file = 0;
  num_receiverCells_total = 0;

  DEALLOCATE_1DPTR(num_receiverCells);
  DEALLOCATE_1DPTR(iblock_receiver);
  DEALLOCATE_1DPTR(ijk_receiverCells);
  DEALLOCATE_1DPTR(iblock_donor);
  DEALLOCATE_1DPTR(ijk_donorCells);
  DEALLOCATE_1DPTR(num_cells_hole);
  DEALLOCATE_1DPTR(ijk_holeCells);

  if (myinput->overset_format == "OVERTURE_NOTSCALING") {

    num_donorBlocks_4myblock = 0;
    num_receiverBlocks_from_myblock = 0;

    DEALLOCATE_1DPTR(iblock_donor_4myblock);
    DEALLOCATE_1DPTR(iblock_receiver_from_myblock);

  } // myinput->overset_format
  else if (myinput->overset_format == "OVERTURE") {

    num_donorRanks_4mygrid = 0;
    num_receiverRanks_from_mygrid = 0;
    total_num_donorCells_4mygrid = 0;
    total_num_receiverCells_from_mygrid = 0;

    DEALLOCATE_1DPTR(irank_receiver);
    DEALLOCATE_1DPTR(irank_donor);

    DEALLOCATE_1DPTR(irank_donor_4mygrid);
    DEALLOCATE_1DPTR(irank_receiver_from_mygrid);

    DEALLOCATE_1DPTR(num_donorCells_4mygrid);
    DEALLOCATE_1DPTR(num_receiverCells_from_mygrid);

  } // myinput->overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // cleanup_overset_temporarydata

} // overture
