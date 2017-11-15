#ifndef IO_IO_H
#define IO_IO_H

//

#include <string>

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../parallel/parallel.h"
#include "../core/input.h"
#include "../geometry/geometry.h"
#include "../core/state.h"
#include "plot3d.h"
#include "overture.h"

namespace io {

class ReportParallel {

  public:

    std::string path;

    void initial_core(void);
    void initial_workLoad(Geometry::StructuredGrid *, int);

}; // ReportParallel



class ReportGeometry {

  public:

    std::string path;

    void initial_region(Geometry::Generic *);
    void initial_block_structured(Geometry::StructuredBlock *);
    void initial_grid_structured(Geometry::StructuredGrid *);

}; // ReportGeometry



class ReportUserInput {

  public:

    std::string path;

    void write_userInput(UserInput);

}; // ReportUserInput



class ReportInitialData {

  public:

    std::string path;

}; // ReportInitialData



class ReportSolution {

  public:

    std::string path;

    void minmax(int, double, Geometry::StructuredGrid *, State *);

}; // ReportSolution



class ReportTiming {

  public:

    std::string path;

}; // ReportTiming



// directory names (make sure to add a slash at the end)
const std::string path_report_todos = "./report/";
const std::string path_report_parallel = "./report/parallel/";
const std::string path_report_geometry = "./report/geometry/";
const std::string path_report_userInput = "./report/userInput/";
const std::string path_report_initialData = "./report/initialData/";
const std::string path_report_solution = "./report/solution/";
const std::string path_report_timing = "./report/timing/";

// file names
const std::string filename_todos = "todos";
const std::string filename_parallel_core = "core";
const std::string filename_parallel_workLoad = "workLoad.dat";
const std::string filename_geometry = "geometry";
const std::string filename_userInput = "userInput";
const std::string filename_initialData = "initialData";
const std::string filename_solution = "solution";
const std::string filename_timing = "timing";

// class objects for reporting
extern ReportParallel report_parallel;
extern ReportGeometry report_geometry;
extern ReportUserInput report_userInput;
extern ReportInitialData report_initialData;
extern ReportSolution report_solution;
extern ReportSolution report_timing;

extern int counter_solution_files;

// basic operations
void create_directory_basic(void);
void mkdir(std::string);
int skip_this_line_of_textFile(std::string);

// standard output
void report_timeadvancing(int, double, double, double, UserInput *, Geometry::StructuredGrid *, State *);

// input deck
const std::string filename_inputDeck = std::string(file_base) + ".in";
void read_inputDeck(int, char * []);
int skip_this_line_of_inputDeck(std::string);

// grid
void read_grid_size(UserInput *, int &, int *&);
void read_grid(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void write_grid(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);

// solution
void read_solution(std::string, UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *);
void read_meanState(std::string, UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, int, double **);
void read_auxvar(std::string, UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, int, double **);
void write_solution(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *);
void write_solution_mean(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *);
void write_solution_aux(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *);
void write_solution_on_the_fly(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *, State *, int);

// metrics
void write_metrics(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);

// overset
void read_overset(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
void cleanup_overset_temporarydata(UserInput *);

// boundary (boundary condition + buffer zone)
struct t_BoundaryData { // useful to contain temporary data read from the boundary condition file

  int iblock; // block ID to which this boundary condition belongs
  int itype; // what kind of boundary treatment
  int idir; // direction of this boundary: XI, ETA, or ZETA
  int iend; // which side of boundary: LEFT or RIGHT
  int is[DIM_MAX], ie[DIM_MAX];

}; // t_BoundaryData
extern t_BoundaryData *boundaryData;
//
void read_bc(UserInput *, Geometry::StructuredGrid *, Geometry::StructuredBlock *);
int skip_this_line_of_bc_file(std::string);
void extract_bcData_out_of_(std::string, Geometry::StructuredGrid *, Geometry::StructuredBlock *, t_BoundaryData &);
int parse_boundary_type(std::string);
int check_if_mygrid_has_this_boundary(Geometry::StructuredGrid *, t_BoundaryData);

// boundary data from separate calculations
extern int num_samples_extern;
extern double period_samples_extern;
extern double *time_extern;
extern double **sol_extern;
//
void read_inflow(UserInput *, Geometry::StructuredGrid *);

void read_file_decompDomain(std::string, int, int []);

} // io

//

#endif
