#ifndef CORE_INPUT_H
#define CORE_INPUT_H

//

#include <string>
#include <vector>

#include "param.h"
#include "macros_inlines.h"
#include "../parallel/parallel.h"

class UserInput {

  public:
    // geometry
    int num_dim;
    int axisym;
    int axis_of_sym;

    // physical model
    std::string model_pde;
    std::string model_fluid;
    std::string model_SGS;

    // passive scalar
    int num_scalar;

    // simulation
    std::string simulation;

    // thermodynamics
    double gamma_specificheat; // ratio of specific heats, gamma = C_p / C_v

    // solution
    int num_vars_sol;
    int num_vars_mean;
    int num_vars_meanGradient;
    int num_vars_aux;

    // temporal discretization
    std::string fix_dt_or_cfl;
    double dt;
    double cfl;
    std::string temporal_scheme;
    int num_time_steps;
    int report_freq;

    // spatial discretization
    std::string spatial_scheme;
    int OA_spatial;
    //
    int do_filter;
    std::string filter_scheme;
    int OA_filter;
    double filter_strength;
    double filter_blend;

    // buffer zone
    int buffer_polynomial_order;
    double buffer_constant;

    // metrics
    std::string scheme_metrics;

    // overset
    int do_overset;
    std::string overset_format;
    std::string overset_type_of_interpolation;
    int overset_accuracy;

    // ghost cells
    int num_cells_ghost;
    int num_cells_ghost_finiteDifference;
    int num_cells_ghost_filter;
    int num_cells_ghost_overset;

    // domain decomposition
    std::string how2decompDomain;
    std::string file_decompDomain;

    // files
    std::string type_file;
    std::string file_grid;
    std::string file_metrics;
    std::string file_solution;
    std::string file_varname_metrics;
    std::string file_varname_solution;

    int present_file_grid_in;
    std::string file_grid_in;

    std::string file_overset;

    int present_file_solution_in;
    std::string file_solution_in;

    int present_file_mean_in;
    std::string file_mean_in;

    int present_file_aux_in;
    std::string file_aux_in;

    std::string file_boundary;

    // solution writing
    double time_writing_solutions;

    // data probing
    int do_probe;
    int num_probes;
    std::string *tmp_probe_name;
    int *tmp_probe_interval;
    double **tmp_probe_xyz;

    // solution interpolation
    std::string interp_fromWhichFormat;
    int num_zones_interpSource;
    int num_dim_interpSource;
    int num_vars_interpSource;
    std::string vars_interpSource;
    double scale_xyz;
    double scale_rho, scale_p, scale_u, scale_s;
    //
    int num_filters_interpolation;
    std::string interpolate_into_which_PLOT3D;
    //
    std::string file_profile_FLUENT;
    std::string file_tecplot_ASCII;

    // time-harmonic wave parameters, if used
    std::string harmonicWave_waveType;
    std::string harmonicWave_waveForm;
    int harmonicWave_idir_propagation;
    double harmonicWave_amplitude;
    double harmonicWave_wavelength, harmonicWave_period;
    double harmonicWave_halfWidth;

    UserInput();
    ~UserInput();

    void set(int, char *[]);
    void set_inputDeck(int, char *[]);
    void set_manual(int, char *[]);

    void check_consistency_dimension(void);
    void check_consistency_between_physical_model_and_simulation(void);
    void check_consistency_wave(void);
    void check_consistency_domainDecomposition(void);

    void get_number_of_ghostCells_due_finiteDifference(void);
    void get_number_of_ghostCells_due_filter(void);
    void get_number_of_ghostCells_due_overset(void);

    void get_number_of_variables(void);

    int truefalse_2int(std::string);
    std::string int_2truefalse(int);
    //
    int xyz_2int(std::string);
    std::string int_2xyz(int);

}; // UserInput

struct t_EntryInputDeck {

  std::string name;
  int id;
  int number_items;
  std::vector<std::string> body;

}; // t_TokenInputDeck

namespace inputDeck {

extern const char delimiter[];

extern int numLinesInputDeck;
extern std::vector<std::string> linesInputDeck;

extern int numEntriesInputDeck;
extern t_EntryInputDeck *entriesInputDeck;

void parse_linesInputDeck(void);
void clear_inputDeck(void);

int count_inputDeck_name(std::string);

int check_inputDeck_name(std::string, int);
int check_inputDeck_keyword(std::string, std::string, int);

void get_userInput(std::string, int &, int = 1);
void get_userInput(std::string, int, int *&, int = 1);
void get_userInput(std::string, double &, int = 1);
void get_userInput(std::string, int, double *&, int = 1);
void get_userInput(std::string, std::string &, int = 1);
void get_userInput(std::string, int, std::vector<std::string> &, int = 1);
//
void get_userInput(std::string, std::string, int &, int = 1);
void get_userInput(std::string, std::string, int, int *&, int = 1);
void get_userInput(std::string, std::string, double &, int = 1);
void get_userInput(std::string, std::string, int, double *&, int = 1);
void get_userInput(std::string, std::string, std::string &, int = 1);
void get_userInput(std::string, std::string, int, std::vector<std::string> &, int = 1);

} // inputDeck

//

#endif
