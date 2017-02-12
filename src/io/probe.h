#ifndef IO_PROBE_H
#define IO_PROBE_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"
#include "../core/input.h"
#include "../geometry/geometry.h"
#include "../core/state.h"

namespace probe {

class ProbePoint {

  public:
    std::string name;
    int interval;
    double xyz[DIM_MAX];
}; // ProbePoint

extern int num_myprobes;
extern ProbePoint *probe_point;

void initialize(UserInput *, Geometry::StructuredGrid *);

} // probe

//

#endif
