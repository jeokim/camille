#ifndef SOLVER_SOURCE_H
#define SOLVER_SOURCE_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"
#include "../geometry/geometry.h"
#include "../core/state.h"

namespace source {

// computes source terms for solution variables
void add_RHSSourceTerms(UserInput *, Geometry::StructuredGrid *, State *, double **, double **, double);

void apply_physicalSource(UserInput *, Geometry::StructuredGrid *, State *, double **, double **, double);

void apply_bufferZone(UserInput *, Geometry::StructuredGrid *, State *, double **, double **, double);
void apply_bufferZone_Freund_ambient(Geometry::StructuredGrid *, State *, double **, double **, double, Geometry::StructuredBufferZone *);
void apply_bufferZone_Freund_Dirichlet(Geometry::StructuredGrid *, State *, double **, double **, double, Geometry::StructuredBufferZone *);
void apply_bufferZone_Freund_harmonicWave(UserInput *, Geometry::StructuredGrid *, State *, double **, double **, double, Geometry::StructuredBufferZone *);
void apply_bufferZone_Freund_file(UserInput *, Geometry::StructuredGrid *, State *, double **, double **, double, Geometry::StructuredBufferZone *);

} // source

//

#endif
