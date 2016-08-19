#ifndef IO_FLUENT_H
#define IO_FLUENT_H

//

#include <string>

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"

namespace fluent {

void read_profile(std::string, int, int, int &, double ** &, double ** &);

} // fluent

//

#endif
