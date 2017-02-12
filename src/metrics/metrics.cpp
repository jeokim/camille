#include "metrics.h"

namespace metrics {

double *x_xi[DIM_MAX][DIM_MAX]; // dx_i/dxi_j (i = 1, 2, 3 for x, y, z; j = 1, 2, 3 for xi, eta, zeta)



void compute(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // initialize
  for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
    for (int irow = XI; irow < DIM_MAX; irow++)
      for (int jcol = XDIR; jcol < DIM_MAX; jcol++)
        mygrid->cell[l0].metrics[irow][jcol] = DUMMY_DOUBLE;

  for (int irow = XDIR; irow < DIM_MAX; irow++)
    for (int jcol = XI; jcol < DIM_MAX; jcol++)
      x_xi[irow][jcol] = new double[mygrid->num_ocells];

  // compute
  if (myinput->scheme_metrics == "THOMAS_LOMBARD") // Thomas & Lombard (AIAA J, 1978)
    compute_ThomasLombard(myinput, mygrid);

  else
    mpi::graceful_exit("Unknown scheme for metrics computation.");

  // Jacobian and inverse Jacobian
  compute_Jacobian(myinput, mygrid);

  // normalize using Jacobian
  normalize_metrics(myinput, mygrid);

  // compute and keep the inverse of normalized metrics
  compute_metricsInverse(myinput, mygrid);

  // clean up
  for (int irow = XDIR; irow < DIM_MAX; irow++)
    for (int jcol = XI; jcol < DIM_MAX; jcol++)
      delete[] x_xi[irow][jcol];

  return;

} // compute



void compute_ThomasLombard(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  int num_cells = mygrid->num_ocells;
  double *func = new double[num_cells];

  for (int idir_phys = XDIR; idir_phys < myinput->num_dim; idir_phys++) {

    // take a derivative in computational coordinates (i.e. xi, eta, and zet)
    for (int idir_comp = XI; idir_comp < myinput->num_dim; idir_comp++) {

      // get a physical coordinate to be differentiated (i.e. x, y, and z)
      for (int l0 = 0; l0 < num_cells; l0++)
        func[l0] = mygrid->cell[l0].xyz[idir_phys];
      spatial::take_derivative_xi_eta_zeta(func, mygrid, idir_comp);

      for (int l0 = 0; l0 < num_cells; l0++)
        (x_xi[idir_phys][idir_comp])[l0] = spatial::dflux[l0];

      // the line below is just to check if dx_i/dxi_j is correctly computed and should not be decommented at production
      //for (int l0 = 0; l0 < num_cells; l0++) mygrid->cell[l0].metrics[idir_phys][idir_comp] = x_xi[idir_phys][idir_comp][l0];

    } // idir_comp
  } // idir_phys

  // dx_i/dxi_j tensor is now available; combine its elements to compute metrics

  // zero out all metrics
  for (int l0 = 0; l0 < num_cells; l0++)
    for (int irow = XI; irow < DIM_MAX; irow++)
      for (int jcol = XDIR; jcol < DIM_MAX; jcol++)
        mygrid->cell[l0].metrics[irow][jcol] = 0.0;

  switch ( myinput->num_dim ) {
  case 1:
    mpi::graceful_exit("1-D metric computation is not implemented yet.");

    break;

  case 2:
    compute_ThomasLombard_2D(myinput, mygrid);

    break;

  case 3:
    compute_ThomasLombard_3D(myinput, mygrid);

    break;

  } // myinput->num_dim

  // clean up
  DEALLOCATE_1DPTR(func);

  return;

} // compute_ThomasLombard



void compute_ThomasLombard_2D(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  int num_cells = mygrid->num_ocells;

  // xi_x
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][XDIR] = (x_xi[YDIR][ETA])[l0];

  // xi_y
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][YDIR] = -(x_xi[XDIR][ETA])[l0];

  // eta_x
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][XDIR] = -(x_xi[YDIR][XI])[l0];

  // eta_y
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][YDIR] = (x_xi[XDIR][XI])[l0];

  return;

} // compute_ThomasLombard_2D



void compute_ThomasLombard_3D(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  int num_cells = mygrid->num_ocells;
  double *func = new double[num_cells];

  // xi_x and zeta_x
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[YDIR][ETA])[l0] * mygrid->cell[l0].xyz[ZDIR]; // y_eta * z
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ZETA); // (y_eta * z)_zeta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][XDIR] = spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, XI); // (y_eta * z)_xi
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ZETA][XDIR] = -spatial::dflux[l0];

  // xi_y and zeta_y
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[XDIR][ETA])[l0] * mygrid->cell[l0].xyz[ZDIR]; // x_eta * z
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ZETA); // (x_eta * z)_zeta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][YDIR] = -spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, XI); // (x_eta * z)_xi
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ZETA][YDIR] = spatial::dflux[l0];

  // xi_z and zeta_z
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[XDIR][ETA])[l0] * mygrid->cell[l0].xyz[YDIR]; // x_eta * y
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ZETA); // (x_eta * y)_zeta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][ZDIR] = spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, XI); // (x_eta * y)_xi
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ZETA][ZDIR] = -spatial::dflux[l0];

  // eta_x and zeta_x
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[YDIR][XI])[l0] * mygrid->cell[l0].xyz[ZDIR]; // y_xi * z
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ZETA); // (y_xi * z)_zeta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][XDIR] = -spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ETA); // (y_xi * z)_eta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ZETA][XDIR] += spatial::dflux[l0];

  // eta_y and zeta_y
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[XDIR][XI])[l0] * mygrid->cell[l0].xyz[ZDIR]; // x_xi * z
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ZETA); // (x_xi * z)_zeta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][YDIR] = spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ETA); // (x_xi * z)_eta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ZETA][YDIR] -= spatial::dflux[l0];

  // eta_z and zeta_z
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[XDIR][XI])[l0] * mygrid->cell[l0].xyz[YDIR]; // x_xi * y
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ZETA); // (x_xi * y)_zeta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][ZDIR] = -spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ETA); // (x_xi * y)_eta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ZETA][ZDIR] += spatial::dflux[l0];

  // xi_x and eta_x
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[YDIR][ZETA])[l0] * mygrid->cell[l0].xyz[ZDIR]; // y_zeta * z
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ETA); // (y_zeta * z)_eta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][XDIR] -= spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, XI); // (y_zeta * z)_xi
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][XDIR] += spatial::dflux[l0];

  // xi_y and eta_y
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[XDIR][ZETA])[l0] * mygrid->cell[l0].xyz[ZDIR]; // x_zeta * z
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ETA); // (x_zeta * z)_eta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][YDIR] += spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, XI); // (x_zeta * z)_xi
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][YDIR] -= spatial::dflux[l0];

  // xi_z and eta_z
  for (int l0 = 0; l0 < num_cells; l0++)
    func[l0] = (x_xi[XDIR][ZETA])[l0] * mygrid->cell[l0].xyz[YDIR]; // x_zeta * y
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, ETA); // (x_zeta * y)_eta
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[XI][ZDIR] -= spatial::dflux[l0];
  //
  spatial::take_derivative_xi_eta_zeta(func, mygrid, XI); // (x_zeta * y)_xi
  for (int l0 = 0; l0 < num_cells; l0++)
    mygrid->cell[l0].metrics[ETA][ZDIR] += spatial::dflux[l0];

  // clean up
  DEALLOCATE_1DPTR(func);

  return;

} // compute_ThomasLombard_3D



void compute_Jacobian(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  //J = det-of-d(xi, eta, zeta)/d(x, y, z)
  //J = (1/J)^{-1} = 1 / det-of-d(x, y, z)/d(xi, eta, zeta)

  int num_cells = mygrid->num_ocells;

  for (int l0 = 0; l0 < num_cells; l0++) {

    mygrid->cell[l0].Jac = DUMMY_DOUBLE;
    mygrid->cell[l0].invJac = DUMMY_DOUBLE;

  } // l0

  // inverse Jacobian first
  switch ( myinput->num_dim ) {
  case 1:

    mpi::graceful_exit("1-D Jacobian computation is not implemented yet.");

    break;

  case 2:
    for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++)
      for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++)
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {

          int l0 = mygrid->idx1D(i, j, k);

          mygrid->cell[l0].invJac = x_xi[XDIR][XI][l0] * x_xi[YDIR][ETA][l0] - x_xi[XDIR][ETA][l0] * x_xi[YDIR][XI][l0];

    } // i

    break;

  case 3:

    mpi::graceful_exit("3-D Jacobian computation is not implemented yet.");
    // following is a straightforward way of computing Jacobian in 3-D
    // however, this suffers from the freestream preservation issue
    // need to be modified to strictly follow Thomas & Lombard (1978)

    //for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++)
    //  for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++)
    //    for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
    //
    //      int l0 = mygrid->idx1D(i, j, k);
    //
    //      mygrid->cell[l0].invJac = x_xi[XDIR][XI][l0] * (x_xi[YDIR][ETA][l0] * x_xi[ZDIR][ZETA][l0] - x_xi[YDIR][ZETA][l0] * x_xi[ZDIR][ETA][l0])
    //                              - x_xi[XDIR][ETA][l0] * (x_xi[YDIR][XI][l0] * x_xi[ZDIR][ZETA][l0] - x_xi[YDIR][ZETA][l0] * x_xi[ZDIR][XI][l0])
    //                              + x_xi[XDIR][ZETA][l0] * (x_xi[YDIR][XI][l0] * x_xi[ZDIR][ETA][l0] - x_xi[YDIR][ETA][l0] * x_xi[ZDIR][XI][l0]);
    //
    //} // i

    break;

  } // myinput->num_dim

  // J = (1/J)^{-1}
  for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++)
    for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++)
      for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {

        int l0 = mygrid->idx1D(i, j, k);

        //assert ( mygrid->cell[l0].invJac > 0.0 );

        mygrid->cell[l0].Jac = 1.0 / mygrid->cell[l0].invJac;

  } // i

  return;

} // compute_Jacobian



void normalize_metrics(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // this function is bypassed since metric is not multiplied with Jacobian in the first place
  // thus, no need to divide by Jacobian here either

  //for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++)
  //  for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++)
  //    for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
  //
  //      int l0 = mygrid->idx1D(i, j, k);
  //
  //      for (int irow = XI; irow < myinput->num_dim; irow++)
  //        for (int jcol = XDIR; jcol < myinput->num_dim; jcol++)
  //          mygrid->cell[l0].metrics[irow][jcol] *= mygrid->cell[l0].invJac;
  //
  //    } // i

  return;

} // normalize_metrics



void compute_metricsInverse(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  // simply invert mygrid->cell[l0].metrics[][]; i.e. (1/J) (dxi_i / dxj)
  // since it is already normalized by Jacobian, no need to attempt to scale the 
  // inverted metrics further

  // simply re-using x_xi[][] is not recommended since metrics may not be simply 
  // its inverse (in 2-D, it could be re-used though); thus, the computed 
  // metrics are inverted for safety

  int num_cells = mygrid->num_ocells;

  // zero out all inverse metrics
  for (int l0 = 0; l0 < num_cells; l0++)
    for (int irow = XDIR; irow < DIM_MAX; irow++)
      for (int jcol = XI; jcol < DIM_MAX; jcol++)
        mygrid->cell[l0].metricsInverse[irow][jcol] = 0.0;

  switch ( myinput->num_dim ) {
  case 1:

    mpi::graceful_exit("1-D metric computation is not implemented yet.");

    break;

  case 2:

    for (int l0 = 0; l0 < num_cells; l0++) {

      mygrid->cell[l0].metricsInverse[XDIR][XI]  =  mygrid->cell[l0].metrics[ETA][YDIR];
      mygrid->cell[l0].metricsInverse[XDIR][ETA] = -mygrid->cell[l0].metrics[XI][YDIR];
      mygrid->cell[l0].metricsInverse[YDIR][XI]  = -mygrid->cell[l0].metrics[ETA][XDIR];
      mygrid->cell[l0].metricsInverse[YDIR][ETA] =  mygrid->cell[l0].metrics[XI][XDIR];

    } // l0

    break;

  case 3:

    for (int l0 = 0; l0 < num_cells; l0++) {

      mygrid->cell[l0].metricsInverse[XDIR][XI]   = mygrid->cell[l0].metrics[ETA][YDIR] * mygrid->cell[l0].metrics[ZETA][ZDIR]
                                                  - mygrid->cell[l0].metrics[ETA][ZDIR] * mygrid->cell[l0].metrics[ZETA][YDIR];
      mygrid->cell[l0].metricsInverse[XDIR][ETA]  = mygrid->cell[l0].metrics[XI][ZDIR]  * mygrid->cell[l0].metrics[ZETA][YDIR]
                                                  - mygrid->cell[l0].metrics[XI][YDIR]  * mygrid->cell[l0].metrics[ZETA][ZDIR];
      mygrid->cell[l0].metricsInverse[XDIR][ZETA] = mygrid->cell[l0].metrics[XI][YDIR]  * mygrid->cell[l0].metrics[ETA][ZDIR]
                                                  - mygrid->cell[l0].metrics[XI][ZDIR]  * mygrid->cell[l0].metrics[ETA][YDIR];

      mygrid->cell[l0].metricsInverse[YDIR][XI]   = mygrid->cell[l0].metrics[ETA][ZDIR] * mygrid->cell[l0].metrics[ZETA][XDIR]
                                                  - mygrid->cell[l0].metrics[ETA][XDIR] * mygrid->cell[l0].metrics[ZETA][ZDIR];
      mygrid->cell[l0].metricsInverse[YDIR][ETA]  = mygrid->cell[l0].metrics[XI][XDIR]  * mygrid->cell[l0].metrics[ZETA][ZDIR]
                                                  - mygrid->cell[l0].metrics[XI][ZDIR]  * mygrid->cell[l0].metrics[ZETA][XDIR];
      mygrid->cell[l0].metricsInverse[YDIR][ZETA] = mygrid->cell[l0].metrics[XI][ZDIR]  * mygrid->cell[l0].metrics[ETA][XDIR]
                                                  - mygrid->cell[l0].metrics[XI][XDIR]  * mygrid->cell[l0].metrics[ETA][ZDIR];

      mygrid->cell[l0].metricsInverse[ZDIR][XI]   = mygrid->cell[l0].metrics[ETA][XDIR] * mygrid->cell[l0].metrics[ZETA][YDIR]
                                                  - mygrid->cell[l0].metrics[ETA][YDIR] * mygrid->cell[l0].metrics[ZETA][XDIR];
      mygrid->cell[l0].metricsInverse[ZDIR][ETA]  = mygrid->cell[l0].metrics[XI][YDIR]  * mygrid->cell[l0].metrics[ZETA][XDIR]
                                                  - mygrid->cell[l0].metrics[XI][XDIR]  * mygrid->cell[l0].metrics[ZETA][YDIR];
      mygrid->cell[l0].metricsInverse[ZDIR][ZETA] = mygrid->cell[l0].metrics[XI][XDIR]  * mygrid->cell[l0].metrics[ETA][YDIR]
                                                  - mygrid->cell[l0].metrics[XI][YDIR]  * mygrid->cell[l0].metrics[ETA][XDIR];

    } // l0

    break;

  } // myinput->num_dim

  return;

} // compute_metricsInverse

} // metrics
