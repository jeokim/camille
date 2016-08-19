#ifndef IO_OVERTURE_H
#define IO_OVERTURE_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../core/input.h"
#include "../parallel/parallel.h"
#include "../geometry/geometry.h"

namespace overture {

extern int num_blocks_in_file;

extern int num_receiverCells_total;
extern int *num_receiverCells;

extern int *iblock_receiver;
extern int *ijk_receiverCells;
//
extern int *iblock_donor;
extern int *ijk_donorCells;
//
extern int *irank_receiver;
extern int *irank_donor;

extern int num_donorBlocks_4myblock;
extern int *iblock_donor_4myblock;
//
extern int num_receiverBlocks_from_myblock;
extern int *iblock_receiver_from_myblock;

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

extern int accuracy_interpolation;
extern int width_stencil[DIM_MAX];
struct t_InterpStencil {

  double *coeff;

}; // t_InterpStencil
extern t_InterpStencil *interpStencilExplicit;

// holes
extern int *num_cells_hole;
extern int *ijk_holeCells;



void read_overset_NOTSCALING(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void read_overset_at_headrank_NOTSCALING(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void propagate_overset_to_grids_NOTSCALING(UserInput *, Geometry::StructuredGrid *);
//
void identify_my_donor_blocks(UserInput *, Geometry::StructuredGrid *);
void identify_my_receiver_blocks(UserInput *, Geometry::StructuredGrid *);
//
void read_overset(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void read_overset_at_headrank(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void propagate_overset_to_grids(UserInput *, Geometry::StructuredGrid *);
//
void convert_blockLevelInfo2GridLevel(UserInput *, Geometry::StructuredGrid *);
void identify_my_donor_ranks(UserInput *, Geometry::StructuredGrid *);
void identify_my_receiver_ranks(UserInput *, Geometry::StructuredGrid *);

void cleanup_overset_temporarydata(UserInput *);

} // overture

//

#endif
