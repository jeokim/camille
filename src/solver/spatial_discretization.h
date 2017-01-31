#ifndef SOLVER_SPATIAL_DISCRETIZATION_H
#define SOLVER_SPATIAL_DISCRETIZATION_H

//

#include <string>

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../core/input.h"
#include "../parallel/parallel.h"
#include "../geometry/geometry.h"
#include "../core/state.h"
#include "../overset/overset.h"
#include "source.h"
#include "boundary_condition.h"

namespace spatial {

extern int num_dim;
extern int num_vars;
extern int num_samples;

extern std::string spatial_scheme;

extern int do_filter; // whether a filter is applied or not
extern std::string filter_scheme; // the type of filtering

extern double *flux;
extern double *dflux;



class SpatialDiscretization {

  public:

    std::string name_scheme;
    int order_of_accuracy;

    // each spatial discretization operator has its own extent over which the operator is applied in my core
    int is_4stencil[DIM_MAX], ie_4stencil[DIM_MAX];
    int iso_4stencil[DIM_MAX], ieo_4stencil[DIM_MAX];

}; // SpatialDiscretization



class StandardCentral : public SpatialDiscretization {

  public:

    int symmetry; // coefficients for boundary stencils are prescribed only on the left end
                  // thus for the right end, coefficients may have an opposite sign depending
                  // on which finite difference I am applying
                  // if symmetry = SYMMETRIC, right-boundary stencil coefficients have the same sign as the left boundary (e.g. filter)
                  // if symmetry = ANTISYMMETRIC, they have the opposite sign (e.g. first derivative)
    int multi_fac_4boundaryRight;

    int half_size_of_stencil;
    double *stencil; // stencil coefficients

    int index_abs(int index_rel) {

      return index_rel + half_size_of_stencil;

    } // index_abs

    int num_of_boundaryCells; // number of cells requiring non-interior stencils
    int size_of_boundaryStencil;
    double *stencil_boundary[MAX_ORDER_ACCURACY / 2]; // for boundary stencils; 
                                                      // N-th order scheme has N/2 boundary points requiring biased or lower-order centered schemes
                                                      // each boundary point will be allocated stencil of which size is size_of_boundaryStencil

    void initialize(UserInput *, Geometry::StructuredGrid *);

    void apply_spatial_discretization(double *, Geometry::StructuredGrid *, int);
    double get_multiplicative_factor_4rightBoundary();

    void take_derivative_xi_eta_zeta(double *, Geometry::StructuredGrid *, int);
    void take_derivative_xyz(double *, Geometry::StructuredGrid *, int);

    // correcting derivatives based upon iblank information
    void iblank_correction(double *, double *, Geometry::StructuredGrid *, int);
    void insert_right_boundary_stencils(Geometry::StructuredGrid *, int *, double *, double *, int);
    void insert_left_boundary_stencils(Geometry::StructuredGrid *, int *, double *, double *, int);
    void check_boundary_stencils_if_all_fluid_points(Geometry::StructuredGrid *, int *, int *, int, int);

}; // StandardCentral



class StandardCentralFilter : public StandardCentral {

  public:

    double filter_strength; // if 0, all-pass filter (minimum strength)
                            // if 1, the highest frequency is completely suppressed (maximum strength)
    double filter_blend; // if 0, filtered variables are not blended with original, unfiltered variables; equivalent to unfiltered case
                         // if 1, filtered variables become the next time-step variables

    void initialize(UserInput *, Geometry::StructuredGrid *);

    void apply_filter(double *, Geometry::StructuredGrid *, int);

}; // StandardCentralFilter



class OptimizedCentralExplicit : public StandardCentral {

  public:

    int num_stencils_interior; // optimized scheme typically doe not have the same order of accuracy as standard scheme of the same stencil size
                               // thus, number of stencil points matter

    void initialize(UserInput *, Geometry::StructuredGrid *);

}; // OptimizedCentralExplicit



class OptimizedCentralFilterExplicit : public StandardCentralFilter {

  public:

    int num_stencils_interior; // optimized scheme typically doe not have the same order of accuracy as standard scheme of the same stencil size
                               // thus, number of stencil points matter

    void initialize(UserInput *, Geometry::StructuredGrid *);

}; // OptimizedCentralFilter



extern StandardCentral cds; // Centered Difference, Standard
extern StandardCentralFilter cfs; // Centered Filter, Standard
extern OptimizedCentralExplicit cdo; // Centered Difference, Optimized
extern OptimizedCentralFilterExplicit cfo; // Centered Filter, Optimized



void initialize(UserInput *myinput, Geometry::StructuredGrid *, State *);
void finalize(void);

int get_size_of_boundaryStencil(void);

// wrapper for taking a derivative: branches depending on which spatial scheme is used
void take_derivative_xi_eta_zeta(double *, Geometry::StructuredGrid *, int);
void take_derivative_xyz(double *, Geometry::StructuredGrid *, int);

// ghost cell operations
void update_ghostcell_data(double *, Geometry::StructuredGrid *, int);
void update_ghostcell_data(double **, int, Geometry::StructuredGrid *, int);
void exchange_ghostcell_data(double *, double *, double *, double *, int, Geometry::StructuredGrid *, int);

// right-hand-side terms of your PDE models; i.e. dy/dt = RHS
void zero_out_flux(void);
void compute_RHS(UserInput *, Geometry::StructuredGrid *, State *, double **, double **, double);
void compute_RHS_acoustics(UserInput *, Geometry::StructuredGrid *, State *, double **, double **);
void compute_RHS_linearized_Euler(UserInput *, Geometry::StructuredGrid *, State *, double **, double **);
void compute_RHS_linearized_Euler_scalar(UserInput *, Geometry::StructuredGrid *, State *, double **, double **);

// update boundary data: boundary condition and overset-grid interpolation
void update_boundary(UserInput *, Geometry::StructuredGrid *, State *, int, double **, double);

// filter
void apply_filter(int, int, double **, Geometry::StructuredGrid *);

void precompute_something(UserInput *, Geometry::StructuredGrid *, State *);

} // spatial

//

#endif
