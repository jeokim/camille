#ifndef CORE_STATE_H
#define CORE_STATE_H

//

#include "param.h"
#include "macros_inlines.h"
#include "input.h"
#include "../math/constants.h"
#include "../parallel/parallel.h"
#include "../geometry/geometry.h"

class State {

  public:
    int num_dim;

    std::string model_pde;

    std::string simulation;

    int num_samples;
    int num_cells_ghost;

    int time_step;
    int time_step_lastrun;
    double time_sol;
    double time_sol_lastrun;

    int num_vars_sol;
    int num_vars_mean;
    int num_vars_meanGradient;
    int num_vars_aux;

    std::string *name_vars;
    std::string *name_vars_mean;
    std::string *name_vars_aux;

    double **sol; // solution vector
    double **sol_old; // solution vector at the last time step
    double *sol_ref; // reference state of solution vector

    double **sol_mean; // mean state of solution vectors
    double **sol_meanGradient; // gradients of mean solution vectors
    double **sol_aux; // auxiliary variables (case dependent)

    double gamma_specificheat;

    void initialize_state(UserInput *, Geometry::StructuredGrid *);

    void backup_the_current_solution(void);

    void compute_dependent_variables(double **);

    void prescribe_on_boundary_solution(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, double, double **);

    // define a one-dimensional mapping of variable indices; useful to store gradients, for example (see sol_meanGradient)
    int ivar1D(int ivar, int idir) {

      return ivar * num_dim + idir;

    } // idx1D

    // compute contravariant velocities
    void to_contravariant_velocity(double (&metrics)[DIM_MAX][DIM_MAX], double velocity_Cartesian[], double velocity_contravariant[]) {

      for (int icontra = XI; icontra < num_dim; icontra++) {

        velocity_contravariant[icontra] = 0.0;
        for (int iCart = XDIR; iCart < num_dim; iCart++)
          velocity_contravariant[icontra] += metrics[icontra][iCart] * velocity_Cartesian[iCart];

      } // icontra

      return;

    } // to_contravariant_velocity

    // compute Cartesian velocities
    void to_Cartesian_velocity(double (&metrics_inverse)[DIM_MAX][DIM_MAX], double Jacobian, double velocity_contravariant[], double velocity_Cartesian[]) {

      for (int iCart = XDIR; iCart < num_dim; iCart++) {

        velocity_Cartesian[iCart] = 0.0;
        for (int icontra = XI; icontra < num_dim; icontra++)
          velocity_Cartesian[iCart] += metrics_inverse[iCart][icontra] * velocity_contravariant[icontra];
        velocity_Cartesian[iCart] *= Jacobian;

      } // icontra

      return;

    } // to_Cartesian_velocity

  private:
    // linear acoustics
    void compute_dependent_variables_acoustics(double **);
    void prescribe_on_boundary_solution_acoustics(Geometry::StructuredBoundaryCondition *, Geometry::StructuredGrid *, double, double **);
    //
    void initialize_state_acoustics(UserInput *, Geometry::StructuredGrid *);
    void initialize_state_acoustics_harmonicWave(UserInput *, Geometry::StructuredGrid *);
    double *acoustics_plane_wave(double, double, int);
    //
    void initialize_state_acoustics_GaussianPulse(UserInput *, Geometry::StructuredGrid *);

    // linearized Euler
    void compute_auxiliary_variables_linearizedEuler(double **);
    void compute_auxiliary_variables_linearizedEuler_mixfrac_constgamma(double **);
    //
    void initialize_state_linearizedEuler(UserInput *, Geometry::StructuredGrid *);
    void initialize_state_linearizedEuler_scalar(UserInput *, Geometry::StructuredGrid *);
    void initialize_state_linearizedEuler_aux_composition(UserInput *, Geometry::StructuredGrid *);
    //
    void initialize_state_linearizedEuler_2Djet_with_a_harmonic_source(UserInput *, Geometry::StructuredGrid *);
    void initialize_state_linearizedEuler_cylinderScattering(UserInput *, Geometry::StructuredGrid *);
    void initialize_state_linearizedEuler_TannaTPN49(UserInput *, Geometry::StructuredGrid *);
    void initialize_state_linearizedEuler_KBKCombustor(UserInput *, Geometry::StructuredGrid *);
    void initialize_state_linearizedEuler_linearNozzle(UserInput *, Geometry::StructuredGrid *);

}; // State

//

#endif
