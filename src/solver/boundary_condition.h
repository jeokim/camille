#ifndef SOLVER_BOUNDARY_CONDITION_H
#define SOLVER_BOUNDARY_CONDITION_H

//

#include "../core/param.h"
#include "../core/input.h"
#include "../parallel/parallel.h"
#include "../geometry/geometry.h"
#include "../core/state.h"
#include "../io/io.h"

namespace bc {

extern int num_vars;
extern int *varIndex_2update;

// following boundary-stencil information is necessary for the Neumann boundary condition
// whose boundary condition depends on boundary finite difference since dq/dxi_i is 
// enforced to the current spatial order of accuracy
extern int num_of_boundaryCells;
extern int size_of_boundaryStencil;
extern double *stencil_boundary[MAX_ORDER_ACCURACY / 2];



void initialize_bc(int, int, int, double *[]);

void apply_BC(UserInput *, Geometry::StructuredGrid *, State *, int, double **, double);

void bc_dirichlet_allzero(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, double, double **);
void bc_dirichlet_harmonicwave(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, State *, double, double **, UserInput *);

void bc_neumann(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, double, double **, int, double **);

void bc_wall_slip_kinematic(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, State *, double, double **, double **);
void bc_wall_slip_kinematic_xyz(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, State *, double, double **, double **, int);

void bc_centerline(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, double, double **, int, int);

} // bc

//

#endif
