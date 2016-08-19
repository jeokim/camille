#include <iostream>
#include <iomanip>

#include "overset.h"

namespace overset {

// overset in general
int do_overset = FALSE;
std::string overset_format = "NONE";
int num_vars_2interpolate = 0;
int accuracy_interpolation = 0;

int num_blocks_in_file;

int num_receiverCells_total;
int *num_receiverCells;
int *ijk_receiverCells;

int *irank_receiver;
int *irank_donor;

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

t_CellReceiver *receiverCells;

int num_donorBlocks_4myblock;
int *iblock_donor_4myblock;

int num_donorCells_in_this_grid = 0;
int width_stencil[DIM_MAX];
t_CellDonor *donorCells;

double *buf_4donors;
int *index4_buf_4donors;

int *ijk_donorCells;

int num_receiverBlocks_from_myblock;
int *iblock_receiver_from_myblock;

t_Buffer4Receivers *buf_4receivers;



void initialize_overset(UserInput *myinput, Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  if (myinput->do_overset == FALSE) return;

  do_overset = myinput->do_overset;
  overset_format = myinput->overset_format;
  num_vars_2interpolate = myinput->num_vars_sol;
  accuracy_interpolation = myinput->overset_accuracy;

  io::read_overset(myinput, mygrid, block);

  num_blocks_in_file = block[0].num_ourkinds;

  initialize_receiver(mygrid, block);
  initialize_donor(mygrid);

  initialize_index4_buf_4donors();

  override_iblank_4holeCells(mygrid);
  override_iblank_4interpolatedCells(mygrid);

  io::cleanup_overset_temporarydata(myinput);
  cleanup_temporarydata();

  mpi::wait_allothers();

  return;

} // initialize_overset



void initialize_receiver(Geometry::StructuredGrid *mygrid, Geometry::StructuredBlock *block) {

  int cur_block = mygrid->id_parent;
  int counter_base;
  int *which_grid, *which_grid_sum;
  int *which_cell, *which_cell_sum;
  int size_of_boundaryStencil = spatial::get_size_of_boundaryStencil();
  int icell;
  //
  int is_block[DIM_MAX], ie_block[DIM_MAX];
  for (int idir = XI; idir < DIM_MAX; idir++) {

    is_block[idir] = block[cur_block].is[idir];
    ie_block[idir] = block[cur_block].ie[idir];

  } // idir

  if (overset_format == "OVERTURE_NOTSCALING") {

    num_receiverCells_total = overture::num_receiverCells_total;

    num_receiverCells = new int[num_blocks_in_file];
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++)
      num_receiverCells[iblock] = overture::num_receiverCells[iblock];

    ijk_receiverCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_receiverCells[counter] = overture::ijk_receiverCells[counter];

    receiverCells = new t_CellReceiver[num_receiverCells[cur_block]];
    for (int counter = 0; counter < num_receiverCells[cur_block]; counter++) {

        receiverCells[counter].igrid = 0;
        receiverCells[counter].l0 = 0;

    } // counter

    counter_base = 0;
    for (int iblock = 0; iblock < cur_block; iblock++)
      counter_base += num_receiverCells[iblock];
    counter_base *= DIM_MAX;
    for (int counter = 0; counter < num_receiverCells[cur_block]; counter++) {

      // get the ijk indices at the block level
      int ijk[DIM_MAX];
      for (int idir = XI; idir < DIM_MAX; idir++)
        ijk[idir] = ijk_receiverCells[counter_base + DIM_MAX*counter + idir];

      // see if my grid contains this cell
      int in_mygrid = FALSE;
      for (int idir = XI; idir < mygrid->num_dim; idir++) {
        if ((ijk[idir] - mygrid->is_in_parent[idir]) * (ijk[idir] - mygrid->ie_in_parent[idir]) <= 0 )
          in_mygrid = TRUE;

        else {

          in_mygrid = FALSE;
          break;

        } // ijk[idir]
      } // idir
      if (in_mygrid == TRUE) {

        receiverCells[counter].igrid = mygrid->id_local;
        receiverCells[counter].l0 = mygrid->idx1D(ijk[XI]   - mygrid->iso_in_parent[XI],  
                                                  ijk[ETA]  - mygrid->iso_in_parent[ETA], 
                                                  ijk[ZETA] - mygrid->iso_in_parent[ZETA]);

//        if (mpi::irank == 2)
//          std::cout << counter << ", ijk: " << std::setw(4) << ijk[XI]
//                               << " "       << std::setw(4) << ijk[ETA]
//                               << " "       << std::setw(4) << ijk[ZETA]
//                               << ", l0: "  << std::setw(4) << receiverCells[counter].l0
//                               << ", xyz: " << std::scientific << mygrid->cell[receiverCells[counter].l0].xyz[XDIR] << " "
//                                            << std::scientific << mygrid->cell[receiverCells[counter].l0].xyz[YDIR] << " "
//                                            << std::scientific << mygrid->cell[receiverCells[counter].l0].xyz[ZDIR] << std::endl;

      } // in_mygrid
//      if (mpi::irank == 2)
//        std::cout << "counter: " << counter << ", which grid?: " << receiverCells[counter].igrid << ", which cell?: " << receiverCells[counter].l0 << std::endl;

    } // counter

    which_grid = new int[num_receiverCells[cur_block]];
    which_grid_sum = new int[num_receiverCells[cur_block]];
    which_cell = new int[num_receiverCells[cur_block]];
    which_cell_sum = new int[num_receiverCells[cur_block]];
    for (int counter = 0; counter < num_receiverCells[cur_block]; counter++) {

      which_grid[counter] = receiverCells[counter].igrid;
      which_grid_sum[counter] = 0.0;
      which_cell[counter] = receiverCells[counter].l0;
      which_cell_sum[counter] = 0.0;

    } // counter
    MPI_Allreduce(which_grid, which_grid_sum, num_receiverCells[cur_block], MPI_INT, MPI_SUM, mpi::comm_block);
    MPI_Allreduce(which_cell, which_cell_sum, num_receiverCells[cur_block], MPI_INT, MPI_SUM, mpi::comm_block);
    for (int counter = 0; counter < num_receiverCells[cur_block]; counter++) {

      receiverCells[counter].igrid = which_grid_sum[counter];
      receiverCells[counter].l0 = which_cell_sum[counter];

    } // counter
    DEALLOCATE_1DPTR(which_grid);
    DEALLOCATE_1DPTR(which_grid_sum);
    DEALLOCATE_1DPTR(which_cell);
    DEALLOCATE_1DPTR(which_cell_sum);
    //
//    if (mpi::irank == 2)
//      for (int counter = 0; counter < num_receiverCells[cur_block]; counter++)
//        std::cout << "counter: " << counter << ", which grid?: " << receiverCells[counter].igrid << ", which cell?: " << receiverCells[counter].l0 << std::endl;

    // store the list of the blocks from which I receive (i.e. my donor blocks)
    num_donorBlocks_4myblock = overture::num_donorBlocks_4myblock;
//    std::cout << "Rank: " << mpi::irank << ", num_donorBlocks_4myblock: " << num_donorBlocks_4myblock << std::endl;
    //
    iblock_donor_4myblock = new int[num_donorBlocks_4myblock];
    for (int iblock = 0; iblock < num_donorBlocks_4myblock; iblock++) {

      iblock_donor_4myblock[iblock] = overture::iblock_donor_4myblock[iblock];
//      std::cout << "Rank: " << mpi::irank << ", iblock: " << iblock << ", iblock_donor_4myblock[iblock]: " << iblock_donor_4myblock[iblock] << std::endl;

    } // iblock

    buf_4donors = new double[num_receiverCells[cur_block] * num_vars_2interpolate];

  } // overset_format
  else if (overset_format == "OVERTURE") {

    num_receiverCells_total = overture::num_receiverCells_total;

    ijk_receiverCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_receiverCells[counter] = overture::ijk_receiverCells[counter];

    irank_receiver = new int[num_receiverCells_total];
    for (int counter = 0; counter < num_receiverCells_total; counter++)
      irank_receiver[counter] = overture::irank_receiver[counter];

    num_donorRanks_4mygrid = overture::num_donorRanks_4mygrid;
    irank_donor_4mygrid = new int[num_donorRanks_4mygrid];
    for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
      irank_donor_4mygrid[irank] = overture::irank_donor_4mygrid[irank];

    num_donorCells_4mygrid = new int[num_donorRanks_4mygrid];
    for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
      num_donorCells_4mygrid[irank] = overture::num_donorCells_4mygrid[irank];
    total_num_donorCells_4mygrid = overture::total_num_donorCells_4mygrid;

    receiverCells = new t_CellReceiver[total_num_donorCells_4mygrid];
    for (int counter = 0; counter < total_num_donorCells_4mygrid; counter++) {

        receiverCells[counter].igrid = 0;
        receiverCells[counter].l0 = 0;

    } // counter

    icell = 0;
    for (int counter = 0; counter < num_receiverCells_total; counter++)
      if (irank_receiver[counter] == mpi::irank) { // if this receiver point belongs to me

        // get the ijk indices at the block level
        int ijk[DIM_MAX];
        for (int idir = XI; idir < DIM_MAX; idir++)
          ijk[idir] = ijk_receiverCells[DIM_MAX*counter + idir];

        receiverCells[icell].igrid = mpi::irank;
        receiverCells[icell].l0 = mygrid->idx1D(ijk[XI]   - mygrid->iso_in_parent[XI],  
                                                ijk[ETA]  - mygrid->iso_in_parent[ETA], 
                                                ijk[ZETA] - mygrid->iso_in_parent[ZETA]);
        icell++;

      } // irank_receiver[counter]
    assert( icell == total_num_donorCells_4mygrid );

    buf_4donors = new double[total_num_donorCells_4mygrid * num_vars_2interpolate];
    index4_buf_4donors = new int[total_num_donorCells_4mygrid];

  } // overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // initialize_receiver



void initialize_donor(Geometry::StructuredGrid *mygrid) {

  int cur_block = mygrid->id_parent;
  int stencilSize;

  if (overset_format == "OVERTURE_NOTSCALING") {

    ijk_donorCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_donorCells[counter] = overture::ijk_donorCells[counter];

    // store the list of the blocks that I donate (i.e. my receiver blocks)
    num_receiverBlocks_from_myblock = overture::num_receiverBlocks_from_myblock;
//    std::cout << "Rank: " << mpi::irank << ", num_receiverBlocks_from_myblock: " << num_receiverBlocks_from_myblock << std::endl;
    //
    iblock_receiver_from_myblock = new int[num_receiverBlocks_from_myblock];
    for (int iblock = 0; iblock < num_receiverBlocks_from_myblock; iblock++) {

      iblock_receiver_from_myblock[iblock] = overture::iblock_receiver_from_myblock[iblock];
//      std::cout << "Rank: " << mpi::irank << ", iblock: " << iblock << ", iblock_receiver_from_myblock[iblock]: " << iblock_receiver_from_myblock[iblock] << std::endl;

    } // iblock

    // allocate buffers to be sent to my receiver blocks
    buf_4receivers = new t_Buffer4Receivers[num_receiverBlocks_from_myblock];
    for (int iblock = 0; iblock < num_receiverBlocks_from_myblock; iblock++) {

      int source_block = cur_block;
      int dest_block = iblock_receiver_from_myblock[iblock];

      int num_cells = num_receiverCells[dest_block];
      int num_vars = num_vars_2interpolate;

      buf_4receivers[iblock].source = source_block;
      buf_4receivers[iblock].dest = dest_block;

      buf_4receivers[iblock].num_cells = num_cells;
      buf_4receivers[iblock].num_vars = num_vars;

      buf_4receivers[iblock].data = new double[num_cells * num_vars];

    } // iblock

    for (int idir = XI; idir < DIM_MAX; idir++)
      width_stencil[idir] = overture::width_stencil[idir];

    // loop over the entire donor points to see if any of those have its block index equal to my block index; 
    // if so, store its corresponding receiver-block index, receiver-cell index, and a grid-level 1-D index
    // note that this grid-level 1-D index corresponds to the lower-left-corner of an interpolation stencil (i.e. smallest i, j, k)
    {
      // first, count num_donorCells_in_this_grid
      num_donorCells_in_this_grid = 0;
      int counter = 0;
      for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
        for (int ircv = 0; ircv < num_receiverCells[iblock]; ircv++) {

          if (overture::iblock_donor[counter] == cur_block) { // at least a grid in our block contains this donor point

            // get the ijk indices at the block level
            int ijk[DIM_MAX];
            for (int idir = XI; idir < DIM_MAX; idir++)
              ijk[idir] = ijk_donorCells[DIM_MAX*counter + idir];

            // see if my grid contains this cell
            int in_mygrid = FALSE;
            for (int idir = XI; idir < mygrid->num_dim; idir++) {
              if ((ijk[idir] - mygrid->is_in_parent[idir]) * (ijk[idir] - mygrid->ie_in_parent[idir]) <= 0 )
                in_mygrid = TRUE;

              else {

                in_mygrid = FALSE;
                break;

              } // ijk[idir]
            } // idir
            if (in_mygrid == TRUE)
              num_donorCells_in_this_grid++;

          } // overture::iblock_donor[counter]
          counter++;

        } // ircv
      } // iblock
      //
      // allocate and store indices
      donorCells = new t_CellDonor[num_donorCells_in_this_grid];
      int counter_donorCells = 0;
      stencilSize = math_algebra::product_array_integer(overture::width_stencil, DIM_MAX);
      counter = 0;
      for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
        for (int ircv = 0; ircv < num_receiverCells[iblock]; ircv++) {

          if (overture::iblock_donor[counter] == cur_block) { // at least a grid in our block contains this donor point

            // get the ijk indices at the block level
            int ijk[DIM_MAX];
            for (int idir = XI; idir < DIM_MAX; idir++)
              ijk[idir] = ijk_donorCells[DIM_MAX*counter + idir];

            // see if my grid contains this cell
            int in_mygrid = FALSE;
            for (int idir = XI; idir < mygrid->num_dim; idir++) {
              if ((ijk[idir] - mygrid->is_in_parent[idir]) * (ijk[idir] - mygrid->ie_in_parent[idir]) <= 0 )
                in_mygrid = TRUE;

              else {

                in_mygrid = FALSE;
                break;

              } // ijk[idir]
            } // idir
            if (in_mygrid == TRUE) {

              donorCells[counter_donorCells].iblock_of_receiverCell = iblock;
              donorCells[counter_donorCells].irank_of_receiverCell = NONE;
              donorCells[counter_donorCells].icell_of_receiverCell = ircv;
              donorCells[counter_donorCells].l0 = mygrid->idx1D(ijk[XI]   - mygrid->iso_in_parent[XI],  
                                                                ijk[ETA]  - mygrid->iso_in_parent[ETA], 
                                                                ijk[ZETA] - mygrid->iso_in_parent[ZETA]);
              donorCells[counter_donorCells].num_coeffs = stencilSize;
              donorCells[counter_donorCells].coeff = new double[stencilSize];
              double sum_coeff = 0.0;
              for (int istencil = 0; istencil < stencilSize; istencil++) {

                donorCells[counter_donorCells].coeff[istencil] = overture::interpStencilExplicit[counter].coeff[istencil];
                sum_coeff += donorCells[counter_donorCells].coeff[istencil];

              } // istencil
              // it is likely that Overture readily did a normalization; but do it here just in case
              sum_coeff = 1.0 / sum_coeff;
              for (int istencil = 0; istencil < stencilSize; istencil++)
                donorCells[counter_donorCells].coeff[istencil] *= sum_coeff;

              counter_donorCells++;

              // check if my grid supports the interpolation stencil in its entirety
              for (int idir = XI; idir < DIM_MAX; idir++) {

                ijk[idir] -= mygrid->iso_in_parent[idir]; // get the local grid-level index of the lower-left-corner of a stencil

                assert( ijk[idir] >= mygrid->iso[idir] ); // make sure the lower-left-corner lies within the left bound of this grid
                assert( ijk[idir] + width_stencil[idir] - 1 <= mygrid->ieo[idir] ); // make sure the right end of the stencil lies within the right bound of this grid
                                                                                    // if this assertion fails, try increasing myinput->num_cells_ghost to accommodate stencils

              } // idir

            } // in_mygrid
          } // overture::iblock_donor[counter]
          counter++;

        } // ircv
      } // iblock
      assert( counter == overture::num_receiverCells_total );
      assert( counter_donorCells == num_donorCells_in_this_grid );

    }

  } // overset_format
  else if (overset_format == "OVERTURE") {

    ijk_donorCells = new int[num_receiverCells_total * DIM_MAX];
    for (int counter = 0; counter < num_receiverCells_total * DIM_MAX; counter++)
      ijk_donorCells[counter] = overture::ijk_donorCells[counter];

    irank_donor = new int[num_receiverCells_total];
    for (int counter = 0; counter < num_receiverCells_total; counter++)
      irank_donor[counter] = overture::irank_donor[counter];

    num_receiverRanks_from_mygrid = overture::num_receiverRanks_from_mygrid;
    irank_receiver_from_mygrid = new int[num_receiverRanks_from_mygrid];
    for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
      irank_receiver_from_mygrid[irank] = overture::irank_receiver_from_mygrid[irank];

    num_receiverCells_from_mygrid = new int[num_receiverRanks_from_mygrid];
    for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
      num_receiverCells_from_mygrid[irank] = overture::num_receiverCells_from_mygrid[irank];
    total_num_receiverCells_from_mygrid = overture::total_num_receiverCells_from_mygrid;

    // allocate buffers to be sent to my receiver blocks
    buf_4receivers = new t_Buffer4Receivers[num_receiverRanks_from_mygrid];
    for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++) {

      int source_rank = mpi::irank;
      int dest_rank = irank_receiver_from_mygrid[irank];

      int num_cells = num_receiverCells_from_mygrid[irank];
      int num_vars = num_vars_2interpolate;

      buf_4receivers[irank].source = source_rank;
      buf_4receivers[irank].dest = dest_rank;

      buf_4receivers[irank].num_cells = num_cells;
      buf_4receivers[irank].num_vars = num_vars;

      buf_4receivers[irank].data = new double[num_cells * num_vars];

    } // irank

    for (int idir = XI; idir < DIM_MAX; idir++)
      width_stencil[idir] = overture::width_stencil[idir];

    // loop over the donor points to see if any of those belong to my grid; 
    // if so, store its corresponding receiver-rank index, receiver-cell index, and a grid-level 1-D index
    // note that this grid-level 1-D index corresponds to the lower-left-corner of an interpolation stencil (i.e. smallest i, j, k)
    {
      // first, count the number of donor cells within this grid
      num_donorCells_in_this_grid = total_num_receiverCells_from_mygrid;

      // allocate and store indices
      donorCells = new t_CellDonor[num_donorCells_in_this_grid];
      int counter_donorCells = 0;
      int counter_receiverCells = 0;
      stencilSize = math_algebra::product_array_integer(overture::width_stencil, DIM_MAX);

      for (int counter = 0; counter < num_receiverCells_total; counter++)
        if (irank_donor[counter] == mpi::irank) { // if this donor point belongs to me

          // get the ijk indices at the block level
          int ijk[DIM_MAX];
          for (int idir = XI; idir < DIM_MAX; idir++)
            ijk[idir] = ijk_donorCells[DIM_MAX*counter + idir];

          // find out the index of irank_receiver_from_mygrid corresponding to 
          // the rank of the receiver point
          int this_index = NONE;
          for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
            if (irank_receiver_from_mygrid[irank] == irank_receiver[counter])
              this_index = irank;
          assert( this_index != NONE );

          donorCells[counter_donorCells].iblock_of_receiverCell = NONE;
          donorCells[counter_donorCells].irank_of_receiverCell = irank_receiver[counter];
          donorCells[counter_donorCells].icell_of_receiverCell = counter_receiverCells;
          donorCells[counter_donorCells].l0 = mygrid->idx1D(ijk[XI]   - mygrid->iso_in_parent[XI],  
                                                            ijk[ETA]  - mygrid->iso_in_parent[ETA], 
                                                            ijk[ZETA] - mygrid->iso_in_parent[ZETA]);
          donorCells[counter_donorCells].num_coeffs = stencilSize;
          donorCells[counter_donorCells].coeff = new double[stencilSize];
          double sum_coeff = 0.0;
          for (int istencil = 0; istencil < stencilSize; istencil++) {

            donorCells[counter_donorCells].coeff[istencil] = overture::interpStencilExplicit[counter].coeff[istencil];
            sum_coeff += donorCells[counter_donorCells].coeff[istencil];

          } // istencil
          // it is likely that Overture readily did a normalization; but do it here just in case
          sum_coeff = 1.0 / sum_coeff;
          for (int istencil = 0; istencil < stencilSize; istencil++)
            donorCells[counter_donorCells].coeff[istencil] *= sum_coeff;

          counter_donorCells++;
          counter_receiverCells++;

          // check if my grid supports the interpolation stencil in its entirety
          for (int idir = XI; idir < DIM_MAX; idir++) {

            ijk[idir] -= mygrid->iso_in_parent[idir]; // get the local grid-level index of the lower-left-corner of a stencil

            assert( ijk[idir] >= mygrid->iso[idir] ); // make sure the lower-left-corner lies within the left bound of this grid
            assert( ijk[idir] + width_stencil[idir] - 1 <= mygrid->ieo[idir] ); // make sure the right end of the stencil lies within the right bound of this grid
                                                                                // if this assertion fails, try increasing myinput->num_cells_ghost to accommodate stencils

          } // idir

      } // counter
      assert( counter_donorCells == num_donorCells_in_this_grid );
      assert( counter_donorCells == counter_receiverCells );

    }

  } // overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // initialize_donor



void initialize_index4_buf_4donors(void) {

  if (overset_format != "OVERTURE")
    return;

  for (int idonor = 0; idonor < total_num_donorCells_4mygrid; idonor++)
    index4_buf_4donors[idonor] = NONE;

  // counter_receiverCells counts receiving cells within my rank
  int counter_receiverCells = 0;

  // counter_donorCells_perRank[] counts cells for the donor ranks sending to me
  int *counter_donorCells_perRank = new int[num_donorRanks_4mygrid];
  for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
    counter_donorCells_perRank[irank] = 0;

  for (int counter = 0; counter < num_receiverCells_total; counter++) {
    if (irank_receiver[counter] == mpi::irank) { // if this receiver point belongs to me

      int this_rank = NONE;
      for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
        if (irank_donor_4mygrid[irank] == irank_donor[counter])
          this_rank = irank;
      assert( this_rank != NONE );

      int counter_base = 0;
      for (int irank = 0; irank < this_rank; irank++)
        counter_base += num_donorCells_4mygrid[irank];
      //
      index4_buf_4donors[counter_base + counter_donorCells_perRank[this_rank]] = counter_receiverCells;

      counter_receiverCells++;
      counter_donorCells_perRank[this_rank]++;

    } // irank_receiver[counter]
  } // counter

  assert( counter_receiverCells == total_num_donorCells_4mygrid );
  //
  int counter = 0;
  for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
    counter += counter_donorCells_perRank[irank];
  assert( counter == total_num_donorCells_4mygrid );

  for (int idonor = 0; idonor < total_num_donorCells_4mygrid; idonor++)
    if (index4_buf_4donors[idonor] < 0 || index4_buf_4donors[idonor] >= total_num_donorCells_4mygrid)
      std::cout << ">>> Something's wrong in initialize_index4_buf_4donors at rank " << mpi::irank 
                << " since index4_buf_4donors[idonor] is " << index4_buf_4donors[idonor] << std::endl;

  DEALLOCATE_1DPTR(counter_donorCells_perRank);

// the implementation below is incorrect since it uses donorCells[idonor].icell_of_receiverCell, 
// which is a donor index, not receiver index (buf_4donors is an array that a receiver core uses).
//  int *indices;
//  for (int irank_send = 0; irank_send < mpi::nprocs; irank_send++) { // loop over the entire cores
//    if ( mpi::irank == irank_send ) { // this is my turn to send the indices
//      for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++) { // loop over the cores I have to send the indices to
//
//        // for simplicity
//        int irank_receiver_thisTurn = irank_receiver_from_mygrid[irank];
//        int num_receiverCells_thisTurn = num_receiverCells_from_mygrid[irank];
//
//        // allocate an array to be sent
//        indices = new int[num_receiverCells_thisTurn];
//        for (int icell = 0; icell < num_receiverCells_thisTurn; icell++)
//          indices[icell] = NONE;
//
//        // fill in the array with donor-cell indices in buf_4donors
//        int counter_donorCells = 0;
//        for (int idonor = 0; idonor < num_donorCells_in_this_grid; idonor++) {
//          if ( donorCells[idonor].irank_of_receiverCell == irank_receiver_thisTurn ) {
//
//            indices[counter_donorCells] = donorCells[idonor].icell_of_receiverCell;
//            counter_donorCells++;
//
//          } // donorCells[idonor].irank_of_receiverCell
//        } // idonor
//        assert( counter_donorCells == num_receiverCells_thisTurn );
//
//        MPI_Send(indices, num_receiverCells_thisTurn, MPI_INT, irank_receiver_thisTurn, mpi::irank, mpi::comm_region);
//
//        DEALLOCATE_1DPTR(indices);
//
//      } // irank
//    } // mpi::irank
//    else { // all the other ranks receive, if necessary
//
//      // check if I am receiving at this turn
//      int this_rank = NONE;
//      for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
//        if (irank_donor_4mygrid[irank] == irank_send)
//          this_rank = irank;
//
//      // allocate and receive
//      if (this_rank != NONE) {
//
//        // for simplicity
//        int num_donorCells_thisTurn = num_donorCells_4mygrid[this_rank];
//
//        // allocate an array to receive
//        indices = new int[num_donorCells_thisTurn];
//
//        MPI_Recv(indices, num_donorCells_thisTurn, MPI_INT, irank_send, irank_send, mpi::comm_region, mpi::status);
//
//        // fill out index4_buf_4donors
//        int counter_base = 0;
//        for (int irank = 0; irank < this_rank; irank++)
//          counter_base += num_donorCells_4mygrid[irank];
//        //
//        for (int idonor = 0; idonor < num_donorCells_thisTurn; idonor++)
//          index4_buf_4donors[counter_base + idonor] = indices[idonor];
//
//        DEALLOCATE_1DPTR(indices);
//
//      } // this_rank
//    } // mpi::irank
//
//  } // irank_send

  return;

} // initialize_index4_buf_4donors



void override_iblank_4holeCells(Geometry::StructuredGrid *mygrid) {

  int cur_block = mygrid->id_parent;
  int counter;
  //
  int *num_cells_hole;
  int *ijk_holeCells;
  int num_cells_hole_total;

  if (overset_format == "OVERTURE_NOTSCALING" || overset_format == "OVERTURE") {

    num_cells_hole = overture::num_cells_hole;
    ijk_holeCells = overture::ijk_holeCells;
    num_cells_hole_total = math_algebra::sum_array_integer(num_cells_hole, num_blocks_in_file);

    counter = 0; // this counts elements in ijk_holeCells
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) { // loop over the entire blocks
      if (num_cells_hole[iblock] > 0) { // if a block contains any hole cells

        if (cur_block == iblock) { // check if the block is my block
          for (int icell = 0; icell < num_cells_hole[iblock]; icell++) {

            // get the ijk indices of a hole cell at the block level
            int ijk[DIM_MAX];
            for (int idir = XI; idir < DIM_MAX; idir++)
              ijk[idir] = overture::ijk_holeCells[idir + counter * DIM_MAX];

            // see if my grid contains this hole cell
            int in_mygrid = FALSE;
            for (int idir = XI; idir < mygrid->num_dim; idir++) {
              if ((ijk[idir] - mygrid->iso_in_parent[idir]) * (ijk[idir] - mygrid->ieo_in_parent[idir]) <= 0 )
                in_mygrid = TRUE;

              else {

                in_mygrid = FALSE;
                break;

              } // ijk[idir]
            } // idir

            if (in_mygrid == TRUE) {

              for (int idir = XI; idir < DIM_MAX; idir++)
                ijk[idir] -= mygrid->iso_in_parent[idir]; // get the local grid-level index of the hole cell

              int l0 = mygrid->idx1D(ijk[XI], ijk[ETA], ijk[ZETA]);
              mygrid->cell[l0].iblank = BLANKED;

            } // in_mygrid

            counter++;

          } // icell
        } // cur_block
        else
          counter += num_cells_hole[iblock];

      } // num_cells_hole[iblock]
    } // iblock
    assert ( counter == num_cells_hole_total );

  } // overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // override_iblank_4holeCells



void override_iblank_4interpolatedCells(Geometry::StructuredGrid *mygrid) {

  int cur_block = mygrid->id_parent;
  int num_cells;

  if (overset_format == "OVERTURE_NOTSCALING") {

    num_cells = num_receiverCells[cur_block];
    for (int icell = 0; icell < num_cells; icell++) {
      if (receiverCells[icell].igrid == mygrid->id_local) {

        int l0 = receiverCells[icell].l0;
        mygrid->cell[l0].iblank = INTERPOLATED;

//        std::cout << "Rank: " << mpi::irank << ", cur_block: " << cur_block << ", xyz: " << mygrid->cell[l0].xyz[XDIR] << " " << mygrid->cell[l0].xyz[YDIR] << " " << mygrid->cell[l0].xyz[ZDIR] << std::endl;

      } // receiverCells[icell].igrid
    } // icell

  } // overset_format
  else if (overset_format == "OVERTURE") {

    for (int icell = 0; icell < total_num_donorCells_4mygrid; icell++) {

      int l0 = receiverCells[icell].l0;
      mygrid->cell[l0].iblank = INTERPOLATED;

    } // icell

  } // overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // override_iblank_4interpolatedCells



void cleanup_temporarydata() {

  if (overset_format == "OVERTURE_NOTSCALING") {

    DEALLOCATE_1DPTR(ijk_receiverCells);
    DEALLOCATE_1DPTR(ijk_donorCells);

  } // overset_format
  else if (overset_format == "OVERTURE") {

    DEALLOCATE_1DPTR(ijk_receiverCells);
    DEALLOCATE_1DPTR(ijk_donorCells);

    DEALLOCATE_1DPTR(irank_receiver);
    DEALLOCATE_1DPTR(irank_donor);

  } // overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // cleanup_temporarydata



void interpolate(Geometry::StructuredGrid *mygrid, double **y) {

  if (overset_format == "OVERTURE_NOTSCALING") {

    empty_buf_4receivers(); // empty buffers before computing the weighted sum
    compute_weighted_sum(mygrid, y); // compute a weighted sum in my grid
    merge_2headranks(mygrid); // collect the contribution from my grids into head ranks
    exchange_at_headranks(mygrid);
    distribute_2mygrids(mygrid);
    update_each_grid(mygrid, y);

  } // overset_format
  else if (overset_format == "OVERTURE") {

    empty_buf_4receivers(); // empty buffers before computing the weighted sum
    compute_weighted_sum(mygrid, y); // compute a weighted sum in my grid
    exchange_at_eachrank(mygrid);
    update_each_grid(mygrid, y);

  } // overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // interpolate



void empty_buf_4receivers() {

  int num_receiverEntities = NONE;

  if (overset_format == "OVERTURE_NOTSCALING") // a head rank of each block collects interpolated
                                               // solutions from its grids and message passing occurs
                                               // among only head ranks; this does not scale with
                                               // the number of cores!
    num_receiverEntities = num_receiverBlocks_from_myblock;

  else if (overset_format == "OVERTURE") // each core does message passing with its donor or receiver pair
    num_receiverEntities = num_receiverRanks_from_mygrid;

  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  for (int ientity = 0; ientity < num_receiverEntities; ientity++) {

    int num_cells = buf_4receivers[ientity].num_cells;
    int num_vars = buf_4receivers[ientity].num_vars;
    int counter = 0;

    for (int icell = 0; icell < num_cells; icell++)
      for (int ivar = 0; ivar < num_vars; ivar++) {

        buf_4receivers[ientity].data[counter] = 0.0;
        counter++;

      } // ivar
  } // ientity

  return;

} // empty_buf_4receivers



void compute_weighted_sum(Geometry::StructuredGrid *mygrid, double **y) {

  int num_vars = num_vars_2interpolate;
  double *interpolant = new double[num_vars];

  int *counter_receiverCells_perRank;
  if (overset_format == "OVERTURE") {

    counter_receiverCells_perRank = new int[num_receiverRanks_from_mygrid];
    for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
      counter_receiverCells_perRank[irank] = 0;

  } // overset_format

  int num_ocells[DIM_MAX];
  for (int idir = XI; idir < DIM_MAX; idir++)
    num_ocells[idir] = mygrid->num_ocells_dir[idir];

  if (num_donorCells_in_this_grid > 0) { // if my grid contains any donor cells
    for (int idonor = 0; idonor < num_donorCells_in_this_grid; idonor++) {

      int iblock_of_receiverCell = donorCells[idonor].iblock_of_receiverCell;
      int irank_of_receiverCell = donorCells[idonor].irank_of_receiverCell;
      int icell_of_receiverCell = donorCells[idonor].icell_of_receiverCell;
      int l0 = donorCells[idonor].l0; // it is my grid-level 1-D index corresponding to the lower-left-corner of the stencil (i.e. smallest i, j, k)
      int num_coeffs = donorCells[idonor].num_coeffs;
      double *coeff = donorCells[idonor].coeff;

      for (int ivar = 0; ivar < num_vars; ivar++)
        interpolant[ivar] = 0.0;
      //
      int istencil = 0;
	    for (int w3 = 0; w3 < width_stencil[ZETA]; w3++)
		    for (int w2 = 0; w2 < width_stencil[ETA]; w2++)
		      for (int w1 = 0; w1 < width_stencil[XI]; w1++) {

            int l1 = l0 + (w1 + w2 * num_ocells[XI] + w3 * num_ocells[XI] * num_ocells[ETA]); // my grid-level 1-D index for stencil members

            for (int ivar = 0; ivar < num_vars; ivar++)
              interpolant[ivar] += coeff[istencil] * (y[ivar])[l1];

            istencil++;

          } // w1

      // find out which buf_4receivers does this donor-cell data belong to
      int which_buf = NONE;
      int ivar_base;
      if (overset_format == "OVERTURE_NOTSCALING") {

        for (int iblock = 0; iblock < num_receiverBlocks_from_myblock; iblock++)
          if (iblock_of_receiverCell == buf_4receivers[iblock].dest)
            which_buf = iblock;
        assert (which_buf != NONE);

        // plug the interpolated solution into a corresponding buf_4receivers
        ivar_base = icell_of_receiverCell * num_vars;
        for (int ivar = 0; ivar < num_vars; ivar++)
           buf_4receivers[which_buf].data[ivar + ivar_base] = interpolant[ivar];

      } // overset_format
      else if (overset_format == "OVERTURE") {

        for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
          if (irank_of_receiverCell == buf_4receivers[irank].dest)
            which_buf = irank;
        assert( which_buf != NONE );

        // plug the interpolated solution into a corresponding buf_4receivers
        ivar_base = counter_receiverCells_perRank[which_buf] * num_vars;
        for (int ivar = 0; ivar < num_vars; ivar++)
           buf_4receivers[which_buf].data[ivar + ivar_base] = interpolant[ivar];

        counter_receiverCells_perRank[which_buf]++;

      } // overset_format
      else
        mpi::graceful_exit("Unknown type of overset-grid format.");

    } // idonor
  } // num_donorCells_in_this_grid

  DEALLOCATE_1DPTR(interpolant);
  if (overset_format == "OVERTURE")
    DEALLOCATE_1DPTR(counter_receiverCells_perRank);

  return;

} // compute_weighted_sum



void merge_2headranks(Geometry::StructuredGrid *mygrid) {

  for (int iblock = 0; iblock < num_receiverBlocks_from_myblock; iblock++) {

    int num_cells = buf_4receivers[iblock].num_cells;
    int num_vars = buf_4receivers[iblock].num_vars;

    double *buf = new double[num_cells * num_vars];
    for (int counter = 0; counter < num_cells * num_vars; counter++)
      buf[counter] = 0.0;

    MPI_Allreduce(buf_4receivers[iblock].data, buf, num_cells * num_vars, MPI_DOUBLE, MPI_SUM, mpi::comm_block);

    for (int counter = 0; counter < num_cells * num_vars; counter++)
      buf_4receivers[iblock].data[counter] = buf[counter];

    DEALLOCATE_1DPTR(buf);

  } // iblock

  return;

} // merge_2headranks



void exchange_at_headranks(Geometry::StructuredGrid *mygrid) {

  int cur_block = mygrid->id_parent;
  double *buf, *buf_all;

  // temporary storage for incoming data from my donor blocks
  if (mpi::irank_block == 0) {

    int num_cells = num_receiverCells[cur_block];
    int num_vars = num_vars_2interpolate;
    int num_data = num_cells * num_vars;

    buf = new double[num_data];
    buf_all = new double[num_donorBlocks_4myblock * num_data];

  } // mpi::irank_block

  // message passing between head ranks
  if (mpi::irank_block == 0) {
    for (int iblock = 0; iblock < num_blocks_in_file; iblock++) {
      if (cur_block == iblock) { // it is my turn to receive from the other blocks

        for (int iblock_donor = 0; iblock_donor < num_donorBlocks_4myblock; iblock_donor++) { // loop over my donors asking to send data

          int num_cells = num_receiverCells[cur_block];
          int num_vars = num_vars_2interpolate;
          int num_data = num_cells * num_vars;
          int counter_ref = iblock_donor * num_data;
          int source = iblock_donor_4myblock[iblock_donor]; // in the MPI group of head ranks (mpi::comm_head), block index is equal to rank

          MPI_Recv(buf, num_data, MPI_DOUBLE, source, source, mpi::comm_head, mpi::status);

          for (int counter = 0; counter < num_data; counter++)
            buf_all[counter + counter_ref] = buf[counter];

        } // iblock_donor
      } // cur_block
      else { // we are sending data now

        for (int iblock_rcv = 0; iblock_rcv < num_receiverBlocks_from_myblock; iblock_rcv++) {

          int num_cells = buf_4receivers[iblock_rcv].num_cells;
          int num_vars = buf_4receivers[iblock_rcv].num_vars;
          int num_data = num_cells * num_vars;
          int dest = buf_4receivers[iblock_rcv].dest; // in the MPI group of head ranks (mpi::comm_head), block index is equal to rank

          if (dest == iblock)
            MPI_Send(buf_4receivers[iblock_rcv].data, num_data, MPI_DOUBLE, dest, mpi::irank_head, mpi::comm_head);

        } // iblock_rcv
      } // cur_block

      mpi::wait_cores_in_(mpi::comm_head);

    } // iblock
  } // mpi::irank_block

  // merge'em
  if (mpi::irank_block == 0) {

    int num_cells = num_receiverCells[cur_block];
    int num_vars = num_vars_2interpolate;
    int num_data = num_cells * num_vars;

    for (int counter = 0; counter < num_data; counter++)
      buf_4donors[counter] = 0.0;

    for (int iblock = 0; iblock < num_donorBlocks_4myblock; iblock++) {

      int counter_ref = iblock * num_data;
      for (int counter = 0; counter < num_data; counter++)
        buf_4donors[counter] += buf_all[counter + counter_ref];

    } // iblock
  } // mpi::irank_block

  // clean up
  if (mpi::irank_block == 0) {

    DEALLOCATE_1DPTR(buf);
    DEALLOCATE_1DPTR(buf_all);

  } // mpi::irank_block

  return;

} // exchange_at_headranks



void exchange_at_eachrank(Geometry::StructuredGrid *mygrid) {

  // each core, if involved in overset-grid interpolation as donor and receiver, 
  // does message passing here with its receiver and donor pair

  int num_vars = num_vars_2interpolate;

  // for efficient message passing, non-blocking send is used; 
  // an implementation using blocking send (commented out below) is a bit slower

// case 1: the following implementation uses non-blocking send

  MPI_Request *request = new MPI_Request[num_receiverRanks_from_mygrid];
  for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++) { // loop over the cores I have to send the interpolated solutions to

    // for simplicity
    int irank_receiver_thisTurn = irank_receiver_from_mygrid[irank];
    int num_receiverCells_thisTurn = num_receiverCells_from_mygrid[irank];
    int num_data = num_receiverCells_thisTurn * num_vars;

//    std::cout << "Sender rank: "         << std::setw(8) << mpi::irank
//              << ", num_receiverRanks: " << std::setw(8) << num_receiverRanks_from_mygrid
//              << ", irank: "             << std::setw(8) << irank
//              << ", irank_receiver: "    << std::setw(8) << irank_receiver_thisTurn
//              << ", num_receiverCells: " << std::setw(8) << num_receiverCells_thisTurn << std::endl;

//    MPI_Send(buf_4receivers[irank].data, num_data, MPI_DOUBLE, irank_receiver_thisTurn, mpi::irank, mpi::comm_region);
    MPI_Isend(buf_4receivers[irank].data, num_data, MPI_DOUBLE, irank_receiver_thisTurn, mpi::irank, mpi::comm_region, &request[irank]);

  } // irank

  int counter_base = 0;
  for (int irank = 0; irank < num_donorRanks_4mygrid; irank++) { // loop over the cores from which the interpolated solutions are sent

    // for simplicity
    int irank_donor_thisTurn = irank_donor_4mygrid[irank];
    int num_donorCells_thisTurn = num_donorCells_4mygrid[irank];
    int num_data = num_donorCells_thisTurn * num_vars;

    // allocate an array and receive
    double *buf = new double[num_data];
    MPI_Recv(buf, num_data, MPI_DOUBLE, irank_donor_thisTurn, irank_donor_thisTurn, mpi::comm_region, mpi::status);

    // place the received data in buf_4donors at correct positions using index4_buf_4donors
    for (int idonor = 0; idonor < num_donorCells_thisTurn; idonor++) {

      int icell = index4_buf_4donors[counter_base + idonor];
      for (int ivar = 0; ivar < num_vars; ivar++)
        buf_4donors[ivar + num_vars * icell] = buf[ivar + num_vars * idonor];

    } // idonor
    counter_base += num_donorCells_thisTurn;

    DEALLOCATE_1DPTR(buf);

  } // irank

  for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++)
    MPI_Wait(&request[irank], mpi::status);
  DEALLOCATE_1DPTR(request);

// case 2: the following implementation uses blocking send
//         it is a little bit less efficient than the block send
//  //
//  for (int irank_send = 0; irank_send < mpi::nprocs; irank_send++) { // loop over the entire cores
//    if ( mpi::irank == irank_send ) { // this is my turn to send the interpolated solutions
//      for (int irank = 0; irank < num_receiverRanks_from_mygrid; irank++) { // loop over the cores I have to send the interpolated solutions to
//
//        // for simplicity
//        int irank_receiver_thisTurn = irank_receiver_from_mygrid[irank];
//        int num_receiverCells_thisTurn = num_receiverCells_from_mygrid[irank];
//        int num_data = num_receiverCells_thisTurn * num_vars;
//
////        std::cout << "Sender rank: "         << std::setw(8) << mpi::irank
////                  << ", num_receiverRanks: " << std::setw(8) << num_receiverRanks_from_mygrid
////                  << ", irank: "             << std::setw(8) << irank
////                  << ", irank_receiver: "    << std::setw(8) << irank_receiver_thisTurn
////                  << ", num_receiverCells: " << std::setw(8) << num_receiverCells_thisTurn << std::endl;
//
//        MPI_Send(buf_4receivers[irank].data, num_data, MPI_DOUBLE, irank_receiver_thisTurn, mpi::irank, mpi::comm_region);
//
//      } // irank
//    } // mpi::irank
//    else { // all the other ranks receive, if necessary
//
//      // check if I am receiving at this turn
//      int this_rank = NONE;
//      for (int irank = 0; irank < num_donorRanks_4mygrid; irank++)
//        if (irank_donor_4mygrid[irank] == irank_send)
//          this_rank = irank;
//
//      // allocate and receive
//      if (this_rank != NONE) {
//
//        // for simplicity
//        int num_donorCells_thisTurn = num_donorCells_4mygrid[this_rank];
//        int num_data = num_donorCells_thisTurn * num_vars;
//
////        sleep(mpi::irank*0.2);
////        std::cout << "   -> Receiver rank: " << std::setw(8) << mpi::irank
////                  << ", num_donorCells: "    << std::setw(8) << num_donorCells_thisTurn << std::endl;
//
//        // allocate an array and receive
//        double *buf = new double[num_data];
//        MPI_Recv(buf, num_data, MPI_DOUBLE, irank_send, irank_send, mpi::comm_region, mpi::status);
//
//        // place the received data in buf_4donors at correct positions using index4_buf_4donors
//        int counter_base = 0;
//        for (int irank = 0; irank < this_rank; irank++)
//          counter_base += num_donorCells_4mygrid[irank];
//        //
//        for (int idonor = 0; idonor < num_donorCells_thisTurn; idonor++) {
//
//          int icell = index4_buf_4donors[counter_base + idonor];
//          for (int ivar = 0; ivar < num_vars; ivar++)
//            buf_4donors[ivar + num_vars * icell] = buf[ivar + num_vars * idonor];
//
//        } // idonor
//
//        DEALLOCATE_1DPTR(buf);
//
//      } // this_rank
//    } // mpi::irank
//  } // irank_send

  return;

} // exchange_at_eachrank



void distribute_2mygrids(Geometry::StructuredGrid *mygrid) {

  // given buf_4donors, simply broadcast within a group consisting of my grids only

  int cur_block = mygrid->id_parent;
  int num_cells = num_receiverCells[cur_block];
  int num_vars = num_vars_2interpolate;

  MPI_Bcast(buf_4donors, num_cells * num_vars, MPI_DOUBLE, 0, mpi::comm_block);

  return;

} // distribute_2mygrids



void update_each_grid(Geometry::StructuredGrid *mygrid, double **y) {

  int cur_block = mygrid->id_parent;
  int num_cells;
  int num_vars = num_vars_2interpolate;

  if (overset_format == "OVERTURE_NOTSCALING") {

    num_cells = num_receiverCells[cur_block];
    for (int icell = 0; icell < num_cells; icell++) {
      if (receiverCells[icell].igrid == mygrid->id_local) {

        int l0 = receiverCells[icell].l0;
        for (int ivar = 0; ivar < num_vars; ivar++)
          (y[ivar])[l0] = buf_4donors[ivar + num_vars * icell];

      } // receiverCells[icell].igrid
    } // icell

  } // overset_format
  else if (overset_format == "OVERTURE") {

    num_cells = total_num_donorCells_4mygrid;
    for (int icell = 0; icell < num_cells; icell++) {

      int l0 = receiverCells[icell].l0;
      for (int ivar = 0; ivar < num_vars; ivar++)
        (y[ivar])[l0] = buf_4donors[ivar + num_vars * icell];

    } // icell

  } // overset_format
  else
    mpi::graceful_exit("Unknown type of overset-grid format.");

  return;

} // update_each_grid

} // overset
