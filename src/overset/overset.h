#ifndef OVERSET_OVERSET_H
#define OVERSET_OVERSET_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../core/input.h"
#include "../geometry/geometry.h"
#include "../core/state.h"
#include "../solver/spatial_discretization.h"
#include "../io/io.h"
#include "../io/overture.h"

namespace overset {

// overset in general
extern int do_overset;
extern std::string overset_format;
extern int num_vars_2interpolate;
extern int accuracy_interpolation;

extern int num_blocks_in_file;

extern int num_receiverCells_total;
extern int *num_receiverCells; // number of receiver cells (cells which are interpolated) per each block
extern int *ijk_receiverCells; // ijk indices of the receivers

extern int *irank_receiver;
extern int *irank_donor;

extern int num_donorRanks_4mygrid;
extern int *irank_donor_4mygrid;
//
extern int num_receiverRanks_from_mygrid;
extern int *irank_receiver_from_mygrid;

extern int *num_donorCells_4mygrid;
extern int total_num_donorCells_4mygrid;
//
extern int *num_receiverCells_from_mygrid;
extern int total_num_receiverCells_from_mygrid;

struct t_CellReceiver {

  int igrid; // grid index to which each receiver (i.e. interpolation point) belongs within a block
  int l0; // the grid-level 1-D index of this receiver cell

}; // t_CellReceiver
extern t_CellReceiver *receiverCells;

extern int num_donorBlocks_4myblock;
extern int *iblock_donor_4myblock;

extern int num_donorCells_in_this_grid;
extern int width_stencil[DIM_MAX];
struct t_CellDonor {

  int iblock_of_receiverCell; // a block index of the corresponding receiver cell
  int irank_of_receiverCell; // a core rank of the corresponding receiver cell
  int icell_of_receiverCell; // a cell index of the corresponding receiver cell
  int l0; // my grid-level 1-D index (not of receiver cell)
  int num_coeffs; // number of stencil weights
  double *coeff; // weights for interpolation

}; // t_CellDonor
extern t_CellDonor *donorCells;

extern double *buf_4donors; // a buffer containing interpolated solutions from my donor blocks
extern int *index4_buf_4donors; // a buffer containing indices for buf_4donors 
                                // data sent from each donor should know where in buf_4donors 
                                // they are inserted; thus, this array stores the indices

extern int *ijk_donorCells; // ijk indices of the donors

extern int num_receiverBlocks_from_myblock;
extern int *iblock_receiver_from_myblock;

struct t_Buffer4Receivers {

  int source;
  int dest;

  int num_cells;
  int num_vars;

  double *data;

}; // t_Buffer4Receivers
extern t_Buffer4Receivers *buf_4receivers;



void initialize_overset(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void initialize_receiver(Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void initialize_donor(Geometry::StructuredGrid *);
void initialize_index4_buf_4donors(void);
void override_iblank_4holeCells(Geometry::StructuredGrid *); // set iblank = BLANKED at hole points
void override_iblank_4interpolatedCells(Geometry::StructuredGrid *); // set iblank = INTERPOLATED at interpolated points
void cleanup_temporarydata(void);

void interpolate(Geometry::StructuredGrid *, double **);
void empty_buf_4receivers(void);
void compute_weighted_sum(Geometry::StructuredGrid *, double **);
void merge_2headranks(Geometry::StructuredGrid *);
void exchange_at_headranks(Geometry::StructuredGrid *);
void exchange_at_eachrank(Geometry::StructuredGrid *);
void distribute_2mygrids(Geometry::StructuredGrid *);
void update_each_grid(Geometry::StructuredGrid *, double **);

} // overset

//

#endif
