#ifndef IO_PLOT3D_H
#define IO_PLOT3D_H

//

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"
#include "../core/input.h"
#include "../geometry/geometry.h"
#include "../core/state.h"

namespace plot3d {

void read_grid_header(std::string, int &, int *&);
void read_grid_header_fileKeptOpen(std::string, int &, int *&, std::ifstream &);
void read_solution_header(std::string, int &, int *&);
void read_solution_header_fileKeptOpen(std::string, int &, int *&, std::ifstream &);
void stdout_grid_size(std::string, int);

void write_grid_header(std::string, Geometry::StructuredBlock *);
void write_solution_header(std::string, Geometry::StructuredBlock *);
void write_function_header(std::string, Geometry::StructuredBlock *, int);

void pack_blockindex_2send(int *, Geometry::StructuredGrid *);
void pack_gridpoint_2send(double *, Geometry::StructuredGrid *);
void pack_iblank_2send(int *, Geometry::StructuredGrid *);
void pack_solution_2send(double *, Geometry::StructuredGrid *, State *);
void pack_function_2send(double *, Geometry::StructuredGrid *, int, double **);

void send_blockindex_2headrank(Geometry::StructuredGrid *, int *, int *);

void read_grid_serialIO(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void minmax_in_gridFile(int, int *, int, double *&, std::ifstream &);
void read_solution_serialIO(std::string, UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *);
void minmax_in_solutionFile(int, int *, int, double *&, double *&, std::ifstream &);
void read_function_serialIO(std::string, UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *, int, double **);
//
void write_grid_serialIO(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void write_solution_serialIO(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *, std::string);
void write_function_serialIO(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, int, double **, std::string);
//
void write_solution_namefile(std::string, State *);
void write_function_namefile(std::string, int, std::string *);
std::string *name_metrics(void);

} // plot3d

//

#endif
