#ifndef PARALLEL_DECOMP_H
#define PARALLEL_DECOMP_H

//

#include <string>

namespace decomp {

void assign_num_cores_crude(int *, int, int *);

void decompDomain(std::string, std::string, int, int [], int, int []);
void decomp_1dir(int [], int, int []);
void decomp_fromFile(std::string, int, int []);

} // decomp


//

#endif
