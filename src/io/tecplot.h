#ifndef IO_TECPLOT_H
#define IO_TECPLOT_H

//

#include <string>

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"

struct t_TecplotFEZone {

  int zoneID;
  int num_dim;
  int num_vars;
  int num_vertices_per_cv;
  int num_nodes;
  int num_cvs;

  double **xyz;
  double **var;
  int **connectivity;

}; // t_TecplotFEZone

namespace tecplot {

extern const int line_width;

void read_ASCII_header_FE(FILE *, int, int);
int read_ASCII_zoneHeader_FE(FILE *, int, int, int &, int &);
void read_ASCII_POINT(FILE *, int, int, int, int, double ** &, double ** &, int ** &, int);

} // tecplot

//

#endif
