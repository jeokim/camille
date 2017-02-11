#ifndef CORE_SIMULATION_H
#define CORE_SIMULATION_H

//

#include "param.h"
#include "macros_inlines.h"
#include "input.h"
#include "../math/constants.h"
#include "../parallel/parallel.h"
#include "../parallel/decomp.h"
#include "../geometry/geometry.h"
#include "state.h"
#include "../io/io.h"
#include "../solver/solver.h"
#include "../metrics/metrics.h"
#include "../overset/overset.h"
#include "../io/probe.h"

//

namespace simulation {

// input
extern UserInput myinput;

// geometry
extern int num_region;
extern int num_blocks_global;
extern int *num_cells_in;
//
extern Geometry::Generic *region;
extern Geometry::StructuredBlock *block;
extern Geometry::StructuredGrid *grid;

// state
extern State state;

void initialize(int, char * []);
void run(void);
void finalize(void);

} // simulation

//

#endif
