#ifndef SOLVER_SOLVER_H
#define SOLVER_SOLVER_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../core/input.h"
#include "../parallel/parallel.h"
#include "../geometry/geometry.h"
#include "../core/state.h"
#include "../io/io.h"
#include "temporal_discretization.h"
#include "spatial_discretization.h"

namespace solver {

void initialize(UserInput *, Geometry::StructuredGrid *, State *);
void finalize();
void precompute_something(UserInput *, Geometry::StructuredGrid *, State *);
void solve(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *);

} // solver

//

#endif
