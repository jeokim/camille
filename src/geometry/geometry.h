#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../core/input.h"
#include "../math/constants.h"
#include "../math/algebra.h"
#include "../parallel/parallel.h"

namespace Geometry {

class GridPoint {

  public:
    double xyz[DIM_MAX];

    double metrics[DIM_MAX][DIM_MAX]; // (1/J) (dxi_i / dxj)
    double metricsInverse[DIM_MAX][DIM_MAX]; // inverse of above; J (dxi / dxi_j)
    double Jac, invJac;

    int iblank; // if BLANKED or 1, this cell is blanked out (masked cell, hole point); i.e. not updated by the flow solver
                // if NOTBLANKED or 0, this cell is a regular fluid cell updated by the flow solver
                // if INTERPOLATED or -1, this cell is interpolated by a grid in another block

}; // GridPoint



class StructuredBoundary {

  public:
    int num_dim;

    int num_cells;
    int num_cells_dir[DIM_MAX];

    int is[DIM_MAX], ie[DIM_MAX];
    int is_in_parent[DIM_MAX], ie_in_parent[DIM_MAX];

    int nXi, nEta, nZeta; // number of cells in each direction used for idx1D()
    void init_idx1D(int, int, int);
    int idx1D(int iXi, int iEta, int iZeta) {

      return nXi * nEta * iZeta + nXi * iEta + iXi; // Xi varies fastest, followed by eta and then zeta

    } // idx1D

    double *xyz[DIM_MAX];

    void initialize(int, int [], int [], int, int);

}; // StructuredBoundary



class StructuredBoundaryCondition : public StructuredBoundary {

  public:
    int which_dir; // along which direction in computational coordinates is this boundary having its normal: e.g. XI or ETA or ZETA
    int which_end; // is it a left boundary or right; e.g. LEFT (with the lowest index; i.e. is) or RIGHT (with the highest index; i.e. ie)
    int which_boundary; // type of this boundary
    int which_model; // which model for the boundary treatment

    void initialize(int, int [], int [], int, int, int);
    void set_type(int, int, int, int);

}; // StructuredBoundaryCondition



class StructuredBufferZone : public StructuredBoundaryCondition {

  public:

    // sponge parameters
    int buffer_polynomial_order;
    double buffer_constant;
    double *buffer_strength;

    StructuredBoundary *boundary;

    void initialize(int, int [], int [], int, int, int);
    void set_parameters_4bufferZone(UserInput *, int [], int []);

}; // StructuredBoundaryCondition



class Generic {

  public:
    int num_dim;

    int num_ourkinds;
    int num_partitions;
    int num_cores;

    int id_global;
    int id_local;

    int id_parent;

    int num_cells;
    int num_ocells;

}; // Generic



class StructuredBlock : public Generic {

  public:
    // ghost cells
    int num_cells_ghost;

    // direction periodicity
    int periodic[DIM_MAX];

    // local indices
    int num_cells_dir[DIM_MAX]; // number of cells
    int num_ocells_dir[DIM_MAX]; // number of physical plus ghost cells
    //
    int is[DIM_MAX], ie[DIM_MAX]; // cell indices
    int iso[DIM_MAX], ieo[DIM_MAX]; // cell indices counting ghost cells

    int nXi, nEta, nZeta; // number of cells in each direction used for idx1D()
    void init_idx1D(int, int, int);
    int idx1D(int iXi, int iEta, int iZeta) { // returns a one-dimensional index for a point in this block
                                              // the returned 1-D index also counts ghost cells

      return nXi * nEta * iZeta + nXi * iEta + iXi; // Xi varies fastest, followed by eta and then zeta

    } // idx1D

    int num_cores_dir[DIM_MAX];

}; // StructuredBlock



class StructuredGrid : public StructuredBlock {

  public:
    // indices in its parent block
    int is_in_parent[DIM_MAX], ie_in_parent[DIM_MAX]; // cell indices
    int iso_in_parent[DIM_MAX], ieo_in_parent[DIM_MAX]; // cell indices counting ghost cells

    int irank_next[DIM_MAX][2]; // local ranks next to me (in both directions)

    // local indices overridden for computing a derivative
    int is_4derivative[DIM_MAX], ie_4derivative[DIM_MAX];
    int iso_4derivative[DIM_MAX], ieo_4derivative[DIM_MAX];

    // local indices overridden for applying a filter
    int is_4filter[DIM_MAX], ie_4filter[DIM_MAX];
    int iso_4filter[DIM_MAX], ieo_4filter[DIM_MAX];

    GridPoint *cell;

    int num_boundaryCondition_nonperiodic;
    StructuredBoundaryCondition *boundaryCondition;

    int num_bufferZones;
    StructuredBufferZone *bufferZone;

    void hardwire_gridPoint(UserInput *, StructuredBlock *);

    // index overriding
    void override_local_indices(UserInput *);
    void override_local_indices_4derivative(UserInput *);
    void override_local_indices_4filter(UserInput *);

    // boundary treatment
    void initialize_boundaryCondition(int);
    void initialize_bufferZone(int);

    // check if a point belongs to this grid
    int check_if_this_is_my_point(int, double [], int *&);

}; // StructuredGrid



// initialize entities
void init_region(Generic *, UserInput *, int, int *);
void init_block(int, StructuredBlock *, Generic *, UserInput *, int [], int);
void init_grid(StructuredGrid *, int, StructuredBlock *, int *);
void init_iblank(StructuredGrid *);

void gen_index_map(StructuredBlock *, int *);

int check_if_this_is_my_block(int, StructuredBlock *);

void clean_up_before_time_marching(StructuredGrid *);

// xyz locations of points in my block
extern double *xyz_block;
void fill_in_xyz_block(double *, int);

// MPI groups for geometric entities
void group_cores_within_region(void);
void group_cores_within_block(StructuredGrid *, StructuredBlock *);
void label_cores_within_block(StructuredGrid *, StructuredBlock *);
void group_cores_head(StructuredBlock *);

void additionalInit_boundary(UserInput *, StructuredGrid *, StructuredBlock *);
void additionalInit_boundaryConditions(UserInput *, StructuredGrid *, StructuredBlock *);
void additionalInit_bufferZones(UserInput *, StructuredGrid *, StructuredBlock *);

} // Geometry

//

#endif
