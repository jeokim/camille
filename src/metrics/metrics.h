#ifndef METRICS_METRICS_H
#define METRICS_METRICS_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"
#include "../geometry/geometry.h"
#include "../solver/spatial_discretization.h"

namespace metrics {

extern double *x_xi[DIM_MAX][DIM_MAX]; // dx_i/dxi_j (i = 1, 2, 3 for x, y, z; j = 1, 2, 3 for xi, eta, zeta)

void compute(UserInput *, Geometry::StructuredGrid *);

void compute_ThomasLombard(UserInput *, Geometry::StructuredGrid *);
void compute_ThomasLombard_2D(UserInput *, Geometry::StructuredGrid *);
void compute_ThomasLombard_3D(UserInput *, Geometry::StructuredGrid *);

void compute_Jacobian(UserInput *, Geometry::StructuredGrid *);

void normalize_metrics(UserInput *, Geometry::StructuredGrid *);

void compute_metricsInverse(UserInput *, Geometry::StructuredGrid *);

} // metrics

//

#endif
