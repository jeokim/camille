#include <cmath>
#include <iostream>
#include <iomanip>

#include "spatial_discretization.h"

namespace spatial {

int num_dim;
int num_vars;
int num_samples;

std::string spatial_scheme;

int do_filter;
std::string filter_scheme;

double *flux;
double *dflux;

StandardCentral cds;
StandardCentralFilter cfs;
OptimizedCentralExplicit cdo;
OptimizedCentralFilterExplicit cfo;



void StandardCentral::initialize(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  name_scheme = "Standard centered finite difference";
  order_of_accuracy = myinput->OA_spatial;
  if (order_of_accuracy%2 != 0)
    mpi::graceful_exit("For centered finite difference, order of accuracy should be even; change it in your input file.");

  // each core has a different extent over which the current spatial discretization is applied
  for (int idir = XI; idir < DIM_MAX; idir++) {

    is_4stencil[idir] = mygrid->is_4derivative[idir];
    ie_4stencil[idir] = mygrid->ie_4derivative[idir];
    iso_4stencil[idir] = mygrid->iso_4derivative[idir];
    ieo_4stencil[idir] = mygrid->ieo_4derivative[idir];

  } // idir

  symmetry = ANTISYMMETRIC;
  multi_fac_4boundaryRight = get_multiplicative_factor_4rightBoundary();

  half_size_of_stencil = order_of_accuracy / 2;
  stencil = new double[order_of_accuracy + 1];

  num_of_boundaryCells = half_size_of_stencil;
  size_of_boundaryStencil = (order_of_accuracy - 2) + 1; // the maximum boundary-stencil size of a N-th order scheme is (N - 2) + 1
  for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) { // N-th order scheme has N/2 boundary points requiring biased or lower-order centered schemes

    stencil_boundary[iboundary] = new double[size_of_boundaryStencil];

    for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++)
      (stencil_boundary[iboundary])[istencil] = 0.0;

  } // iboundary

  // See Anderson, Tannehill, Pletcher (Taylor & Francis) or Fornberg (Mathematics of Computation 1988)
  switch ( order_of_accuracy ) {
  case 2: // second-order standard centered finite difference

    // interior stencil
    stencil[index_abs(0)] =  0.0;
    stencil[index_abs(1)] =  1.0 / 2.0;

    // i = 0: first-order biased
    stencil_boundary[FIRST][0] = -1.0; // corresponding to i = 0
    stencil_boundary[FIRST][1] =  1.0;

    break;

  case 4: // fourth-order standard centered finite difference

    // interior stencil
    stencil[index_abs(0)] =  0.0;
    stencil[index_abs(1)] =  8.0 / 12.0;
    stencil[index_abs(2)] = -1.0 / 12.0;

    // i = 0: second-order biased
    stencil_boundary[FIRST][0] = -3.0 / 2.0; // corresponding to i = 0
    stencil_boundary[FIRST][1] =  4.0 / 2.0;
    stencil_boundary[FIRST][2] = -1.0 / 2.0;

    // i = 1: second-order centered
    stencil_boundary[SECOND][0] = -1.0 / 2.0;
    stencil_boundary[SECOND][1] =  0.0; // corresponding to i = 1
    stencil_boundary[SECOND][2] =  1.0 / 2.0;

    break;

  default:

    mpi::graceful_exit("Unknown order of accuracy for standard centered finite difference; implement its coefficients first.");

    break;

  } // order_of_accuracy

  // ensure interior-stencil anti-symmetry so that the finite difference is non-dissipative
  for (int istencil = 1; istencil <= half_size_of_stencil; istencil++)
    stencil[index_abs(-istencil)] = -stencil[index_abs(istencil)];

  return;

} // StandardCentral::initialize



void StandardCentralFilter::initialize(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  name_scheme = "Standard centered filter";
  order_of_accuracy = myinput->OA_filter;
  if (order_of_accuracy%2 != 0)
    mpi::graceful_exit("For centered filter, order of accuracy should be even; change it in your input file.");
  if (order_of_accuracy > MAX_ORDER_ACCURACY)
    mpi::graceful_exit("MAX_ORDER_ACCURACY is exceeded for filter.");

  // each core has a different extent over which the current spatial discretization is applied
  for (int idir = XI; idir < DIM_MAX; idir++) {

    is_4stencil[idir] = mygrid->is_4filter[idir];
    ie_4stencil[idir] = mygrid->ie_4filter[idir];
    iso_4stencil[idir] = mygrid->iso_4filter[idir];
    ieo_4stencil[idir] = mygrid->ieo_4filter[idir];

  } // idir

  symmetry = SYMMETRIC;
  multi_fac_4boundaryRight = get_multiplicative_factor_4rightBoundary();

  half_size_of_stencil = order_of_accuracy / 2;
  stencil = new double[order_of_accuracy + 1];

  num_of_boundaryCells = half_size_of_stencil;
  size_of_boundaryStencil = (order_of_accuracy - 2) + 1; // the maximum boundary-stencil size of a N-th order scheme is (N - 2) + 1
  for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) { // N-th order scheme has N/2 boundary points requiring biased or lower-order centered schemes

    stencil_boundary[iboundary] = new double[size_of_boundaryStencil];

    for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++)
      (stencil_boundary[iboundary])[istencil] = 0.0;

  } // iboundary

  filter_strength = myinput->filter_strength; // note that standard centeral filter has always the maximum strength; i.e. 1.0
  filter_blend = myinput->filter_blend;

  switch ( order_of_accuracy ) {
  case 4: // fourth-order standard centered filter

    // note that both interior stencil and boundary cells are filtered

    // interior stencil; fourth-order standard centered filter; see C.2.4 of Lele (JCP 1992) with alpha = 0 = beta and d = 0
    stencil[index_abs(0)] = 5.0 / 8.0;
    stencil[index_abs(1)] = (1.0 / 2.0) / 2.0;
    stencil[index_abs(2)] = (-1.0 / 8.0) / 2.0;

    // i = 0: not filtered
    stencil_boundary[FIRST][0] = 1.0; // corresponding to i = 0

    // i = 1: not filtered
    stencil_boundary[SECOND][1] = 1.0; // corresponding to i = 1

    break;

  case 8: // eighth-order standard centered filter

    // note that interior stencil is defined using damping function and boundary cells are filtered

    // interior stencil; see SFs9p on p212 of Bogey & Bailly (JCP 2004); the filter formula is found in equation (2)
    stencil[index_abs(0)] = 35.0 / 128.0;
    stencil[index_abs(1)] = -7.0 / 32.0;
    stencil[index_abs(2)] =  7.0 / 64.0;
    stencil[index_abs(3)] = -1.0 / 32.0;
    stencil[index_abs(4)] =  1.0 / 256.0;

    // include the filter strength in the coefficients
    stencil[index_abs(0)] = 1.0 - filter_strength * stencil[index_abs(0)];
    for (int istencil = 1; istencil <= half_size_of_stencil; istencil++)
      stencil[index_abs(istencil)] *= -filter_strength;

    // i = 0: not filtered
    stencil_boundary[FIRST][0] = 1.0; // corresponding to i = 0

    // i = 1: not filtered
    stencil_boundary[SECOND][1] = 1.0; // corresponding to i = 1

    // i = 2: fourth-order standard centered filter; see C.2.4 of Lele (JCP 1992) with alpha = 0 = beta and d = 0
    stencil_boundary[THIRD][0] = (-1.0 / 8.0) / 2.0;
    stencil_boundary[THIRD][1] = (1.0 / 2.0) / 2.0;
    stencil_boundary[THIRD][2] = 5.0 / 8.0; // corresponding to i = 2
    stencil_boundary[THIRD][3] = (1.0 / 2.0) / 2.0;
    stencil_boundary[THIRD][4] = (-1.0 / 8.0) / 2.0;

    // i = 3: sixth-order standard centered filter; see C.2.5 of Lele (JCP 1992)
    stencil_boundary[FOURTH][0] = (1.0 / 32.0) / 2.0;
    stencil_boundary[FOURTH][1] = (-3.0 / 16.0) / 2.0;
    stencil_boundary[FOURTH][2] = (15.0 / 32.0) / 2.0;
    stencil_boundary[FOURTH][3] = 11.0 / 16.0; // corresponding to i = 3
    stencil_boundary[FOURTH][4] = (15.0 / 32.0) / 2.0;
    stencil_boundary[FOURTH][5] = (-3.0 / 16.0) / 2.0;
    stencil_boundary[FOURTH][6] = (1.0 / 32.0) / 2.0;

    break;

  default:

    mpi::graceful_exit("Unknown order of accuracy for standard centered filter; implement its coefficients first.");

    break;

  } // order_of_accuracy

  // ensure interior-stencil symmetry so that the filter is purely dissipative
  for (int istencil = 1; istencil <= half_size_of_stencil; istencil++)
    stencil[index_abs(-istencil)] = stencil[index_abs(istencil)];

  return;

} // StandardCentralFilter::initialize



void OptimizedCentralExplicit::initialize(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  name_scheme = "Optimized centered finite difference, explicit";
  order_of_accuracy = myinput->OA_spatial;
  if (order_of_accuracy%2 != 0)
    mpi::graceful_exit("For centered finite difference, order of accuracy should be even; change it in your input file.");

  // each core has a different extent over which the current spatial discretization is applied
  for (int idir = XI; idir < DIM_MAX; idir++) {

    is_4stencil[idir] = mygrid->is_4derivative[idir];
    ie_4stencil[idir] = mygrid->ie_4derivative[idir];
    iso_4stencil[idir] = mygrid->iso_4derivative[idir];
    ieo_4stencil[idir] = mygrid->ieo_4derivative[idir];

  } // idir

  // for optimized scheme, its formal order of accuracy does not describe much about the stencil
  num_stencils_interior = NONE;
  if (myinput->spatial_scheme == "FDO11P")
    num_stencils_interior = 11;
  else
    mpi::graceful_exit("Unknown optimized finite difference scheme.");
  //
  if (num_stencils_interior > MAX_ORDER_ACCURACY + 1)
    mpi::graceful_exit("MAX_ORDER_ACCURACY is exceeded for finite difference.");

  symmetry = ANTISYMMETRIC;
  multi_fac_4boundaryRight = get_multiplicative_factor_4rightBoundary();

  half_size_of_stencil = (num_stencils_interior - 1) / 2;
  stencil = new double[num_stencils_interior];

  num_of_boundaryCells = half_size_of_stencil;
  size_of_boundaryStencil = num_stencils_interior - 2; // e.g. if 11-point scheme, a biggest boundary stencil comes from a 9-point scheme to preserve symmetry even for boundary scheme
  for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) { // N-point scheme has N/2 boundary points requiring biased or lower-order centered schemes

    stencil_boundary[iboundary] = new double[size_of_boundaryStencil];

    for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++)
      (stencil_boundary[iboundary])[istencil] = 0.0;

  } // iboundary

  switch ( num_stencils_interior ) {
  case 11: // eleven-point optimized finite difference, explicit

    // interior stencil; see FDo11p on p212 of Bogey & Bailly (JCP 2004); the finite difference formula is found in equation (1)
    stencil[index_abs(0)] =  0.0;
    stencil[index_abs(1)] =  0.872756993962;
    stencil[index_abs(2)] = -0.286511173973;
    stencil[index_abs(3)] =  0.090320001280;
    stencil[index_abs(4)] = -0.020779405824;
    stencil[index_abs(5)] =  0.002484594688;

    // i = 0: fourth order seven-point FD06 (Berland et al., JCP 2007)
    stencil_boundary[FIRST][0] = -2.225833963270; // corresponding to i = 0
    stencil_boundary[FIRST][1] =  4.827779580575;
    stencil_boundary[FIRST][2] = -5.001388453836;
    stencil_boundary[FIRST][3] =  3.911103941646;
    stencil_boundary[FIRST][4] = -2.115267458633;
    stencil_boundary[FIRST][5] =  0.718882784412;
    stencil_boundary[FIRST][6] = -0.115276430894;

    // i = 1: fourth order seven-point FD15 (Berland et al., JCP 2007)
    stencil_boundary[SECOND][0] = -0.212932721951;
    stencil_boundary[SECOND][1] = -1.060320390770; // corresponding to i = 1
    stencil_boundary[SECOND][2] =  2.078926116439;
    stencil_boundary[SECOND][3] = -1.287179452384;
    stencil_boundary[SECOND][4] =  0.685176395471;
    stencil_boundary[SECOND][5] = -0.245320613994;
    stencil_boundary[SECOND][6] =  0.041650667189;

    // i = 2: fourth order five-point standard scheme
    stencil_boundary[THIRD][0] =  1.0 / 12.0;
    stencil_boundary[THIRD][1] = -8.0 / 12.0;
    stencil_boundary[THIRD][2] =  0.0; // corresponding to i = 2
    stencil_boundary[THIRD][3] =  8.0 / 12.0;
    stencil_boundary[THIRD][4] = -1.0 / 12.0;

    // i = 3: fourth order seven-point optimized FDo7p of Bogey & Bailly (JCP 2004)
    stencil_boundary[FOURTH][0] = -0.02625099337517627;
    stencil_boundary[FOURTH][1] =  0.18833730683403838;
    stencil_boundary[FOURTH][2] = -0.797921633542548;
    stencil_boundary[FOURTH][3] =  0.0; // corresponding to i = 3
    stencil_boundary[FOURTH][4] =  0.797921633542548;
    stencil_boundary[FOURTH][5] = -0.18833730683403838;
    stencil_boundary[FOURTH][6] =  0.02625099337517627;

    // i = 4: fourth order nine-point optimized FDo9p of Bogey & Bailly (JCP 2004)
    stencil_boundary[FIFTH][0] =  0.007650904064;
    stencil_boundary[FIFTH][1] = -0.059463584768;
    stencil_boundary[FIFTH][2] =  0.244678631765;
    stencil_boundary[FIFTH][3] = -0.841570125482;
    stencil_boundary[FIFTH][4] =  0.0; // corresponding to i = 4
    stencil_boundary[FIFTH][5] =  0.841570125482;
    stencil_boundary[FIFTH][6] = -0.244678631765;
    stencil_boundary[FIFTH][7] =  0.059463584768;
    stencil_boundary[FIFTH][8] = -0.007650904064;

    break;

  default:

    mpi::graceful_exit("Unknown number of stencils for optimized centered finite difference, explicit; implement its coefficients first.");

    break;

  } // num_stencils_interior

  // ensure interior-stencil anti-symmetry so that the finite difference is non-dissipative
  for (int istencil = 1; istencil <= half_size_of_stencil; istencil++)
    stencil[index_abs(-istencil)] = -stencil[index_abs(istencil)];

  return;

} // OptimizedCentralExplicit::initialize



void OptimizedCentralFilterExplicit::initialize(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  name_scheme = "Optimized centered filter, explicit";
  order_of_accuracy = myinput->OA_filter;
  if (order_of_accuracy%2 != 0)
    mpi::graceful_exit("For centered filter, order of accuracy should be even; change it in your input file.");

  // each core has a different extent over which the current spatial discretization is applied
  for (int idir = XI; idir < DIM_MAX; idir++) {

    is_4stencil[idir] = mygrid->is_4filter[idir];
    ie_4stencil[idir] = mygrid->ie_4filter[idir];
    iso_4stencil[idir] = mygrid->iso_4filter[idir];
    ieo_4stencil[idir] = mygrid->ieo_4filter[idir];

  } // idir

  // for optimized scheme, its formal order of accuracy does not describe much about the stencil
  num_stencils_interior = NONE;
  if (myinput->filter_scheme == "SFO11P")
    num_stencils_interior = 11;
  else
    mpi::graceful_exit("Unknown optimized filter scheme.");
  //
  if (num_stencils_interior > MAX_ORDER_ACCURACY + 1)
    mpi::graceful_exit("MAX_ORDER_ACCURACY is exceeded for filter.");

  symmetry = SYMMETRIC;
  multi_fac_4boundaryRight = get_multiplicative_factor_4rightBoundary();

  half_size_of_stencil = (num_stencils_interior - 1) / 2;
  stencil = new double[num_stencils_interior];

  num_of_boundaryCells = half_size_of_stencil;
  size_of_boundaryStencil = num_stencils_interior - 2; // e.g. if 11-point scheme, a biggest boundary stencil comes from a 9-point scheme to preserve symmetry even for boundary scheme
  for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) { // N-point scheme has N/2 boundary points requiring biased or lower-order centered schemes

    stencil_boundary[iboundary] = new double[size_of_boundaryStencil];

    for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++)
      (stencil_boundary[iboundary])[istencil] = 0.0;

  } // iboundary

  filter_strength = myinput->filter_strength;
  filter_blend = myinput->filter_blend;

  switch ( num_stencils_interior ) {
  case 11: // eleven-point optimized centered filter, explicit

    // note that only interior stencil and i = 4 are defined in terms of damping function

    // interior stencil; see SFo11p on p212 of Bogey & Bailly (JCP 2004); the filter formula is found in equation (2)
    stencil[index_abs(0)] =  0.215044884112;
    stencil[index_abs(1)] = -0.187772883589;
    stencil[index_abs(2)] =  0.123755948787;
    stencil[index_abs(3)] = -0.059227575576;
    stencil[index_abs(4)] =  0.018721609157;
    stencil[index_abs(5)] = -0.002999540835;

    // include the filter strength in the coefficients
    stencil[index_abs(0)] = 1.0 - filter_strength * stencil[index_abs(0)];
    for (int istencil = 1; istencil <= half_size_of_stencil; istencil++)
      stencil[index_abs(istencil)] *= -filter_strength;

    // i = 0: not filtered
    stencil_boundary[FIRST][0] = 1.0; // corresponding to i = 0

    // i = 1: not filtered
    stencil_boundary[SECOND][1] = 1.0; // corresponding to i = 1

    // i = 2: fourth-order standard centered filter; see C.2.4 of Lele (JCP 1992) with alpha = 0 = beta and d = 0
    stencil_boundary[THIRD][0] = (-1.0 / 8.0) / 2.0;
    stencil_boundary[THIRD][1] = (1.0 / 2.0) / 2.0;
    stencil_boundary[THIRD][2] = 5.0 / 8.0; // corresponding to i = 2
    stencil_boundary[THIRD][3] = (1.0 / 2.0) / 2.0;
    stencil_boundary[THIRD][4] = (-1.0 / 8.0) / 2.0;

    // i = 3: sixth-order standard centered filter; see C.2.5 of Lele (JCP 1992)
    stencil_boundary[FOURTH][0] = (1.0 / 32.0) / 2.0;
    stencil_boundary[FOURTH][1] = (-3.0 / 16.0) / 2.0;
    stencil_boundary[FOURTH][2] = (15.0 / 32.0) / 2.0;
    stencil_boundary[FOURTH][3] = 11.0 / 16.0; // corresponding to i = 3
    stencil_boundary[FOURTH][4] = (15.0 / 32.0) / 2.0;
    stencil_boundary[FOURTH][5] = (-3.0 / 16.0) / 2.0;
    stencil_boundary[FOURTH][6] = (1.0 / 32.0) / 2.0;

//    // i = 0: not damped
//    stencil_boundary[FIRST][0] = 0.0; // corresponding to i = 0
//
//    // i = 1: not damped
//    stencil_boundary[SECOND][1] = 0.0; // corresponding to i = 1
//
//    // i = 2: not damped
//    stencil_boundary[THIRD][2] = 0.0; // corresponding to i = 2
//
//    // i = 3: fourth-order seven-point optimized filter, explicit (SFo7p)
//    stencil_boundary[FOURTH][0] = -0.0257391577927332;
//    stencil_boundary[FOURTH][1] =  0.1139783155854664;
//    stencil_boundary[FOURTH][2] = -0.2242608422072667;
//    stencil_boundary[FOURTH][3] =  0.2720433688290671; // corresponding to i = 3
//    stencil_boundary[FOURTH][4] = -0.2242608422072667;
//    stencil_boundary[FOURTH][5] =  0.1139783155854664;
//    stencil_boundary[FOURTH][6] = -0.0257391577927332;
//    // i = 3: not damped
//    stencil_boundary[FOURTH][3] = 0.0; // corresponding to i = 3

    // i = 4: fourth-order nine-point optimized filter, explicit; see SFo9p on p212 of Bogey & Bailly (JCP 2004)
    stencil_boundary[FIFTH][0] =  0.008228661760;
    stencil_boundary[FIFTH][1] = -0.045211119360;
    stencil_boundary[FIFTH][2] =  0.120007591680;
    stencil_boundary[FIFTH][3] = -0.204788880640;
    stencil_boundary[FIFTH][4] =  0.243527493120; // corresponding to i = 4
    stencil_boundary[FIFTH][5] = -0.204788880640;
    stencil_boundary[FIFTH][6] =  0.120007591680;
    stencil_boundary[FIFTH][7] = -0.045211119360;
    stencil_boundary[FIFTH][8] =  0.008228661760;
    //
    for (int iboundary = FIFTH; iboundary <= FIFTH; iboundary++)
      for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++) {

        if (iboundary == istencil)
          stencil_boundary[iboundary][istencil] = 1.0 - filter_strength * stencil_boundary[iboundary][istencil];
        else
          stencil_boundary[iboundary][istencil] *= -filter_strength;

      } // istencil

    break;

  default:

    mpi::graceful_exit("Unknown number of stencils for optimized centered filter, explicit; implement its coefficients first.");

    break;

  } // num_stencils_interior

  // ensure interior-stencil symmetry so that the filter is purely dissipative
  for (int istencil = 1; istencil <= half_size_of_stencil; istencil++)
    stencil[index_abs(-istencil)] = stencil[index_abs(istencil)];

  return;

} // OptimizedCentralFilterExplicit::initialize



int get_size_of_boundaryStencil() {

  int size_of_boundaryStencil;

  if (spatial_scheme == "STANDARD_CENTRAL")
    size_of_boundaryStencil = cds.size_of_boundaryStencil;

  else if (spatial_scheme == "FDO11P")
    size_of_boundaryStencil = cdo.size_of_boundaryStencil;

  else
    mpi::graceful_exit("Unknown spatial scheme.");

  return size_of_boundaryStencil;

} // get_size_of_boundaryStencil



void update_ghostcell_data(double *func, Geometry::StructuredGrid *mygrid, int idir) {

  int source_left = mygrid->irank_next[idir][LEFT];
  int source_right = mygrid->irank_next[idir][RIGHT];
  //
  int np[DIM_MAX-1]; // number of points on a plane where grid indices in idir are equal
  np[FIRST] = mygrid->num_cells_dir[dir_other[idir][FIRST]];
  np[SECOND] = mygrid->num_cells_dir[dir_other[idir][SECOND]];
  //
  int count = mygrid->num_cells_ghost * np[FIRST] * np[SECOND];
  double *dbuf_send_left = new double[count];
  double *dbuf_send_right = new double[count];
  double *dbuf_recv_left = new double[count];
  double *dbuf_recv_right = new double[count];
  for (int i = 0; i < count; i++) {

    dbuf_send_left[i] = DUMMY_DOUBLE;
    dbuf_send_right[i] = DUMMY_DOUBLE;
    dbuf_recv_left[i] = DUMMY_DOUBLE;
    dbuf_recv_right[i] = DUMMY_DOUBLE;

  } // i

  int N1, N12;
  switch ( idir ) {
  case XI:
    // pack donor points in my grid
    N1 = mygrid->num_cells_ghost;
    N12 = N1 * mygrid->num_cells_dir[ETA];
    //
    for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
      int kk = k - mygrid->is[ZETA];
      for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
        int jj = j - mygrid->is[ETA];
        for (int i_stencil = 0; i_stencil < mygrid->num_cells_ghost; i_stencil++) {

          int l1 = kk * N12 + jj * N1 + i_stencil;

          int i = mygrid->is[idir] + i_stencil;
          int l0 = mygrid->idx1D(i, j, k);
          dbuf_send_left[l1] = func[l0]; // data going to the core on the left
          //
          i = mygrid->ie[idir] - mygrid->num_cells_ghost + 1 + i_stencil;
          l0 = mygrid->idx1D(i, j, k);
          dbuf_send_right[l1] = func[l0]; // data going to the core on the right

        } // i_stencil
      } // j
    } // k

    exchange_ghostcell_data(dbuf_send_left, dbuf_send_right, dbuf_recv_left, dbuf_recv_right, count, mygrid, idir);

    // unpack neighbor's if there is any
    if (source_right != NONE) {
      for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
        int kk = k - mygrid->is[ZETA];
        for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
          int jj = j - mygrid->is[ETA];
          for (int i_stencil = 0; i_stencil < mygrid->num_cells_ghost; i_stencil++) {

            int l1 = kk * N12 + jj * N1 + i_stencil;

            int i = mygrid->ie[idir] + 1 + i_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            func[l0] = dbuf_recv_right[l1]; // data coming from the core on the right

          } // i_stencil
        } // j
      } // k
    } // source_right
    if (source_left != NONE) {
      for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
        int kk = k - mygrid->is[ZETA];
        for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
          int jj = j - mygrid->is[ETA];
          for (int i_stencil = 0; i_stencil < mygrid->num_cells_ghost; i_stencil++) {

            int l1 = kk * N12 + jj * N1 + i_stencil;

            int i = mygrid->is[idir] - mygrid->num_cells_ghost + i_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            func[l0] = dbuf_recv_left[l1]; // data coming from the core on the left

          } // i_stencil
        } // j
      } // k
    } // source_left

    break;

  case ETA:
    // pack donor points in my grid
    N1 = mygrid->num_cells_ghost;
    N12 = N1 * mygrid->num_cells_dir[XI];
    //
    for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
      int kk = k - mygrid->is[ZETA];
      for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
        int ii = i - mygrid->is[XI];
        for (int j_stencil = 0; j_stencil < mygrid->num_cells_ghost; j_stencil++) {

          int l1 = kk * N12 + ii * N1 + j_stencil;

          int j = mygrid->is[idir] + j_stencil;
          int l0 = mygrid->idx1D(i, j, k);
          dbuf_send_left[l1] = func[l0]; // data going to the core on the left
          //
          j = mygrid->ie[idir] - mygrid->num_cells_ghost + 1 + j_stencil;
          l0 = mygrid->idx1D(i, j, k);
          dbuf_send_right[l1] = func[l0]; // data going to the core on the right

        } // j_stencil
      } // i
    } // k

    exchange_ghostcell_data(dbuf_send_left, dbuf_send_right, dbuf_recv_left, dbuf_recv_right, count, mygrid, idir);

    // unpack neighbor's if there is any
    if (source_right != NONE) {
      for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
        int kk = k - mygrid->is[ZETA];
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
          int ii = i - mygrid->is[XI];
          for (int j_stencil = 0; j_stencil < mygrid->num_cells_ghost; j_stencil++) {

            int l1 = kk * N12 + ii * N1 + j_stencil;

            int j = mygrid->ie[idir] + 1 + j_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            func[l0] = dbuf_recv_right[l1]; // data coming from the core on the right

          } // j_stencil
        } // i
      } // k
    } // source_right
    if (source_left != NONE) {
      for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
        int kk = k - mygrid->is[ZETA];
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
          int ii = i - mygrid->is[XI];
          for (int j_stencil = 0; j_stencil < mygrid->num_cells_ghost; j_stencil++) {

            int l1 = kk * N12 + ii * N1 + j_stencil;

            int j = mygrid->is[idir] - mygrid->num_cells_ghost + j_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            func[l0] = dbuf_recv_left[l1]; // data coming from the core on the left

          } // j_stencil
        } // i
      } // k
    } // source_left

    break;

  case ZETA:
    // pack donor points in my grid
    N1 = mygrid->num_cells_ghost;
    N12 = N1 * mygrid->num_cells_dir[XI];
    //
    for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
      int jj = j - mygrid->is[ETA];
      for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
        int ii = i - mygrid->is[XI];
        for (int k_stencil = 0; k_stencil < mygrid->num_cells_ghost; k_stencil++) {

          int l1 = jj * N12 + ii * N1 + k_stencil;

          int k = mygrid->is[idir] + k_stencil;
          int l0 = mygrid->idx1D(i, j, k);
          dbuf_send_left[l1] = func[l0]; // data going to the core on the left
          //
          k = mygrid->ie[idir] - mygrid->num_cells_ghost + 1 + k_stencil;
          l0 = mygrid->idx1D(i, j, k);
          dbuf_send_right[l1] = func[l0]; // data going to the core on the right

        } // k_stencil
      } // i
    } // j

    exchange_ghostcell_data(dbuf_send_left, dbuf_send_right, dbuf_recv_left, dbuf_recv_right, count, mygrid, idir);

    // unpack neighbor's if there is any
    if (source_right != NONE) {
      for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
        int jj = j - mygrid->is[ETA];
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
          int ii = i - mygrid->is[XI];
          for (int k_stencil = 0; k_stencil < mygrid->num_cells_ghost; k_stencil++) {

            int l1 = jj * N12 + ii * N1 + k_stencil;

            int k = mygrid->ie[idir] + 1 + k_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            func[l0] = dbuf_recv_right[l1]; // data coming from the core on the right

          } // k_stencil
        } // i
      } // j
    } // source_right
    if (source_left != NONE) {
      for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
        int jj = j - mygrid->is[ETA];
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
          int ii = i - mygrid->is[XI];
          for (int k_stencil = 0; k_stencil < mygrid->num_cells_ghost; k_stencil++) {

            int l1 = jj * N12 + ii * N1 + k_stencil;

            int k = mygrid->is[idir] - mygrid->num_cells_ghost + k_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            func[l0] = dbuf_recv_left[l1]; // data coming from the core on the left

          } // k_stencil
        } // i
      } // j
    } // source_left

    break;

  } // idir

  // clean up
  DEALLOCATE_1DPTR(dbuf_send_left);
  DEALLOCATE_1DPTR(dbuf_send_right);
  DEALLOCATE_1DPTR(dbuf_recv_left);
  DEALLOCATE_1DPTR(dbuf_recv_right);

  return;

} // update_ghostcell_data



void update_ghostcell_data(double **func, int num_vars_in, Geometry::StructuredGrid *mygrid, int idir) {

  // this is a modified version of 
  //   void update_ghostcell_data(double *func, Geometry::StructuredGrid *mygrid, int idir)
  // instead of a one-dimensional aray "double *func", a two-dimensional array "double **func" 
  // with the number of variables, "num_vars_in" are tossed in
  // the way it works is essentially the same as its original one-dimensional case

  int source_left = mygrid->irank_next[idir][LEFT];
  int source_right = mygrid->irank_next[idir][RIGHT];
  //
  int np[DIM_MAX-1]; // number of points on a plane where grid indices in idir are equal
  np[FIRST] = mygrid->num_cells_dir[dir_other[idir][FIRST]];
  np[SECOND] = mygrid->num_cells_dir[dir_other[idir][SECOND]];
  //
  int count = num_vars_in * (mygrid->num_cells_ghost * np[FIRST] * np[SECOND]);
  double *dbuf_send_left = new double[count];
  double *dbuf_send_right = new double[count];
  double *dbuf_recv_left = new double[count];
  double *dbuf_recv_right = new double[count];
  for (int i = 0; i < count; i++) {

    dbuf_send_left[i] = DUMMY_DOUBLE;
    dbuf_send_right[i] = DUMMY_DOUBLE;
    dbuf_recv_left[i] = DUMMY_DOUBLE;
    dbuf_recv_right[i] = DUMMY_DOUBLE;

  } // i

  int N1, N12, N123;
  switch ( idir ) {
  case XI:
    // pack donor points in my grid
    N1 = mygrid->num_cells_ghost;
    N12 = N1 * mygrid->num_cells_dir[ETA];
    N123 = N12 * mygrid->num_cells_dir[ZETA];
    //
    for (int ivar = 0; ivar < num_vars_in; ivar++) {
      for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
        int kk = k - mygrid->is[ZETA];
        for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
          int jj = j - mygrid->is[ETA];
          for (int i_stencil = 0; i_stencil < mygrid->num_cells_ghost; i_stencil++) {

            int l1 = ivar * N123 + kk * N12 + jj * N1 + i_stencil;

            int i = mygrid->is[idir] + i_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            dbuf_send_left[l1] = (func[ivar])[l0]; // data going to the core on the left
            //
            i = mygrid->ie[idir] - mygrid->num_cells_ghost + 1 + i_stencil;
            l0 = mygrid->idx1D(i, j, k);
            dbuf_send_right[l1] = (func[ivar])[l0]; // data going to the core on the right

          } // i_stencil
        } // j
      } // k
    } // ivar

    exchange_ghostcell_data(dbuf_send_left, dbuf_send_right, dbuf_recv_left, dbuf_recv_right, count, mygrid, idir);

    // unpack neighbor's if there is any
    if (source_right != NONE) {
      for (int ivar = 0; ivar < num_vars_in; ivar++) {
        for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
          int kk = k - mygrid->is[ZETA];
          for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
            int jj = j - mygrid->is[ETA];
            for (int i_stencil = 0; i_stencil < mygrid->num_cells_ghost; i_stencil++) {

              int l1 = ivar * N123 + kk * N12 + jj * N1 + i_stencil;

              int i = mygrid->ie[idir] + 1 + i_stencil;
              int l0 = mygrid->idx1D(i, j, k);
              (func[ivar])[l0] = dbuf_recv_right[l1]; // data coming from the core on the right

            } // i_stencil
          } // j
        } // k
      } // ivar
    } // source_right
    if (source_left != NONE) {
      for (int ivar = 0; ivar < num_vars_in; ivar++) {
        for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
          int kk = k - mygrid->is[ZETA];
          for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
            int jj = j - mygrid->is[ETA];
            for (int i_stencil = 0; i_stencil < mygrid->num_cells_ghost; i_stencil++) {

              int l1 = ivar * N123 + kk * N12 + jj * N1 + i_stencil;

              int i = mygrid->is[idir] - mygrid->num_cells_ghost + i_stencil;
              int l0 = mygrid->idx1D(i, j, k);
              (func[ivar])[l0] = dbuf_recv_left[l1]; // data coming from the core on the left

            } // i_stencil
          } // j
        } // k
      } // ivar
    } // source_left

    break;

  case ETA:
    // pack donor points in my grid
    N1 = mygrid->num_cells_ghost;
    N12 = N1 * mygrid->num_cells_dir[XI];
    N123 = N12 * mygrid->num_cells_dir[ZETA];
    //
    for (int ivar = 0; ivar < num_vars_in; ivar++) {
      for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
        int kk = k - mygrid->is[ZETA];
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
          int ii = i - mygrid->is[XI];
          for (int j_stencil = 0; j_stencil < mygrid->num_cells_ghost; j_stencil++) {

            int l1 = ivar * N123 + kk * N12 + ii * N1 + j_stencil;

            int j = mygrid->is[idir] + j_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            dbuf_send_left[l1] = (func[ivar])[l0]; // data going to the core on the left
            //
            j = mygrid->ie[idir] - mygrid->num_cells_ghost + 1 + j_stencil;
            l0 = mygrid->idx1D(i, j, k);
            dbuf_send_right[l1] = (func[ivar])[l0]; // data going to the core on the right

          } // j_stencil
        } // i
      } // k
    } // ivar

    exchange_ghostcell_data(dbuf_send_left, dbuf_send_right, dbuf_recv_left, dbuf_recv_right, count, mygrid, idir);

    // unpack neighbor's if there is any
    if (source_right != NONE) {
      for (int ivar = 0; ivar < num_vars_in; ivar++) {
        for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
          int kk = k - mygrid->is[ZETA];
          for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
            int ii = i - mygrid->is[XI];
            for (int j_stencil = 0; j_stencil < mygrid->num_cells_ghost; j_stencil++) {

              int l1 = ivar * N123 + kk * N12 + ii * N1 + j_stencil;

              int j = mygrid->ie[idir] + 1 + j_stencil;
              int l0 = mygrid->idx1D(i, j, k);
              (func[ivar])[l0] = dbuf_recv_right[l1]; // data coming from the core on the right

            } // j_stencil
          } // i
        } // k
      } // ivar
    } // source_right
    if (source_left != NONE) {
      for (int ivar = 0; ivar < num_vars_in; ivar++) {
        for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++) {
          int kk = k - mygrid->is[ZETA];
          for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
            int ii = i - mygrid->is[XI];
            for (int j_stencil = 0; j_stencil < mygrid->num_cells_ghost; j_stencil++) {

              int l1 = ivar * N123 + kk * N12 + ii * N1 + j_stencil;

              int j = mygrid->is[idir] - mygrid->num_cells_ghost + j_stencil;
              int l0 = mygrid->idx1D(i, j, k);
              (func[ivar])[l0] = dbuf_recv_left[l1]; // data coming from the core on the left

            } // j_stencil
          } // i
        } // k
      } // ivar
    } // source_left

    break;

  case ZETA:
    // pack donor points in my grid
    N1 = mygrid->num_cells_ghost;
    N12 = N1 * mygrid->num_cells_dir[XI];
    N123 = N12 * mygrid->num_cells_dir[ETA];
    //
    for (int ivar = 0; ivar < num_vars_in; ivar++) {
      for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
        int jj = j - mygrid->is[ETA];
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
          int ii = i - mygrid->is[XI];
          for (int k_stencil = 0; k_stencil < mygrid->num_cells_ghost; k_stencil++) {

            int l1 = ivar * N123 + jj * N12 + ii * N1 + k_stencil;

            int k = mygrid->is[idir] + k_stencil;
            int l0 = mygrid->idx1D(i, j, k);
            dbuf_send_left[l1] = (func[ivar])[l0]; // data going to the core on the left
            //
            k = mygrid->ie[idir] - mygrid->num_cells_ghost + 1 + k_stencil;
            l0 = mygrid->idx1D(i, j, k);
            dbuf_send_right[l1] = (func[ivar])[l0]; // data going to the core on the right

          } // k_stencil
        } // i
      } // j
    } // ivar

    exchange_ghostcell_data(dbuf_send_left, dbuf_send_right, dbuf_recv_left, dbuf_recv_right, count, mygrid, idir);

    // unpack neighbor's if there is any
    if (source_right != NONE) {
      for (int ivar = 0; ivar < num_vars_in; ivar++) {
        for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
          int jj = j - mygrid->is[ETA];
          for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
            int ii = i - mygrid->is[XI];
            for (int k_stencil = 0; k_stencil < mygrid->num_cells_ghost; k_stencil++) {

              int l1 = ivar * N123 + jj * N12 + ii * N1 + k_stencil;

              int k = mygrid->ie[idir] + 1 + k_stencil;
              int l0 = mygrid->idx1D(i, j, k);
              (func[ivar])[l0] = dbuf_recv_right[l1]; // data coming from the core on the right

            } // k_stencil
          } // i
        } // j
      } // ivar
    } // source_right
    if (source_left != NONE) {
      for (int ivar = 0; ivar < num_vars_in; ivar++) {
        for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++) {
          int jj = j - mygrid->is[ETA];
          for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
            int ii = i - mygrid->is[XI];
            for (int k_stencil = 0; k_stencil < mygrid->num_cells_ghost; k_stencil++) {

              int l1 = ivar * N123 + jj * N12 + ii * N1 + k_stencil;

              int k = mygrid->is[idir] - mygrid->num_cells_ghost + k_stencil;
              int l0 = mygrid->idx1D(i, j, k);
              (func[ivar])[l0] = dbuf_recv_left[l1]; // data coming from the core on the left

            } // k_stencil
          } // i
        } // j
      } // ivar
    } // source_left

    break;

  } // idir

  // clean up
  DEALLOCATE_1DPTR(dbuf_send_left);
  DEALLOCATE_1DPTR(dbuf_send_right);
  DEALLOCATE_1DPTR(dbuf_recv_left);
  DEALLOCATE_1DPTR(dbuf_recv_right);

  return;

} // update_ghostcell_data



void exchange_ghostcell_data(double *dbuf_send_left, double *dbuf_send_right, double *dbuf_recv_left, double *dbuf_recv_right, 
                             int count, Geometry::StructuredGrid *mygrid, int idir) {

  int source, dest;

  // message passing occurs among cores within a same block
  for (int irank_recv = 0; irank_recv < mpi::nprocs_block; irank_recv++) {
    if (mpi::irank_block == irank_recv) { // if this is my turn to receive

      // receive from the right
      source = mygrid->irank_next[idir][RIGHT];
      if (source != NONE) // receive from the core on my right, if exists
        MPI_Recv(dbuf_recv_right, count, MPI_DOUBLE, source, source, mpi::comm_block, mpi::status);

      // receive from the left
      source = mygrid->irank_next[idir][LEFT];
      if (source != NONE) // receive from the core on my left, if exists
        MPI_Recv(dbuf_recv_left, count, MPI_DOUBLE, source, source, mpi::comm_block, mpi::status);

    } // mpi::irank_block
    else {

      // send to the left
      dest = mygrid->irank_next[idir][LEFT];
      if (dest == irank_recv) // send to the core on my left, if it is a correct receiver
        MPI_Send(dbuf_send_left, count, MPI_DOUBLE, dest, mpi::irank_block, mpi::comm_block);

      // send to the right
      dest = mygrid->irank_next[idir][RIGHT];
      if (dest == irank_recv) // send to the core on my right, if it is a correct receiver
        MPI_Send(dbuf_send_right, count, MPI_DOUBLE, dest, mpi::irank_block, mpi::comm_block);

    } // mpi::irank_block

//    mpi::wait_cores_in_(mpi::comm_block); // do not remove the comment here!

  } // irank_recv

//  // send to the left
//  int dest = mygrid->irank_next[idir][LEFT];
//  if (dest != NONE) // send to the core on my left, if exists
//    MPI_Isend(dbuf_send_left, count, MPI_DOUBLE, dest, mpi::irank_block, mpi::comm_block, &mpi::request[0]);
//  //
//  // receive from the right
//  int source = mygrid->irank_next[idir][RIGHT];
//  if (source != NONE) // receive from the core on my right, if exists
//    MPI_Irecv(dbuf_recv_right, count, MPI_DOUBLE, source, source, mpi::comm_block, &mpi::request[1]);
//  //
//  mpi::wait_cores_in_(mpi::comm_block);
//
//  // send to the right
//  dest = mygrid->irank_next[idir][RIGHT];
//  if (dest != NONE) // send to the core on my right, if exists
//    MPI_Isend(dbuf_send_right, count, MPI_DOUBLE, dest, mpi::irank_block, mpi::comm_block, &mpi::request[0]);
//  //
//  // receive from the left
//  source = mygrid->irank_next[idir][LEFT];
//  if (source != NONE) // receive from the core on my left, if exists
//    MPI_Irecv(dbuf_recv_left, count, MPI_DOUBLE, source, source, mpi::comm_block, &mpi::request[1]);
//  //
//  mpi::wait_cores_in_(mpi::comm_block);

  return;

} // exchange_ghostcell_data



void StandardCentral::apply_spatial_discretization(double *func, Geometry::StructuredGrid *mygrid, int idir_operator) {

  if (mygrid->num_cells_dir[idir_operator] == 1) // OPERATOR_IF_ONECELL
    return;

  // initialization
  for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
    dflux[l0] = 0.0;

  // update ghost cells
  update_ghostcell_data(func, mygrid, idir_operator);

  switch ( idir_operator ) {
  case XI:

    // interior scheme
    for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++)
      for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++)
        for (int i = is_4stencil[XI]; i <= ie_4stencil[XI]; i++) {

          int l0 = mygrid->idx1D(i, j, k);

          for (int l = -half_size_of_stencil; l <= half_size_of_stencil; l++) {

            int l1 = l0 + l;
            dflux[l0] += stencil[index_abs(l)] * func[l1];

          } // l
        } // i

    // boundary scheme
    // left boundary
    if (mygrid->irank_next[idir_operator][LEFT] == NONE) {
      for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++)
        for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++) {

          int lso = mygrid->idx1D(iso_4stencil[XI], j, k);
          for (int i = iso_4stencil[XI]; i < is_4stencil[XI]; i++) {
            int ilocal = i - iso_4stencil[XI];

            int l0 = mygrid->idx1D(i, j, k);

            for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++) {

              int l1 = lso + istencil;
              dflux[l0] += stencil_boundary[ilocal][istencil] * func[l1];

            } // istencil
          } // i
        } // j
    } // mygrid->irank_next[idir_operator][LEFT]

    // right boundary
    if (mygrid->irank_next[idir_operator][RIGHT] == NONE) {
      for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++)
        for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++) {

          int leo = mygrid->idx1D(ieo_4stencil[XI], j, k);
          for (int i = ieo_4stencil[XI]; i > ie_4stencil[XI]; i--) {
            int ilocal = ieo_4stencil[XI] - i;

            int l0 = mygrid->idx1D(i, j, k);

            for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++) {

              int l1 = leo - istencil;
              dflux[l0] += multi_fac_4boundaryRight * stencil_boundary[ilocal][istencil] * func[l1];

            } // istencil
          } // i
        } // j
    } // mygrid->irank_next[idir_operator][RIGHT]

    break;

  case ETA:

    // interior scheme
    for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++)
      for (int j = is_4stencil[ETA]; j <= ie_4stencil[ETA]; j++)
        for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {

          int l0 = mygrid->idx1D(i, j, k);

          for (int l = -half_size_of_stencil; l <= half_size_of_stencil; l++) {

            int l1 = l0 + mygrid->nXi * l;
            dflux[l0] += stencil[index_abs(l)] * func[l1];

          } // l
        } // i

    // boundary scheme
    // left boundary
    if (mygrid->irank_next[idir_operator][LEFT] == NONE) {
      for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++)
        for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {

          int lso = mygrid->idx1D(i, iso_4stencil[ETA], k);
          for (int j = iso_4stencil[ETA]; j < is_4stencil[ETA]; j++) {
            int jlocal = j - iso_4stencil[ETA];

            int l0 = mygrid->idx1D(i, j, k);

            for (int jstencil = 0; jstencil < size_of_boundaryStencil; jstencil++) {

              int l1 = lso + mygrid->nXi * jstencil;
              dflux[l0] += stencil_boundary[jlocal][jstencil] * func[l1];

            } // jstencil
          } // j
        } // i
    } // mygrid->irank_next[idir_operator][LEFT]

    // right boundary
    if (mygrid->irank_next[idir_operator][RIGHT] == NONE) {
      for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++)
        for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {

          int leo = mygrid->idx1D(i, ieo_4stencil[ETA], k);
          for (int j = ieo_4stencil[ETA]; j > ie_4stencil[ETA]; j--) {
            int jlocal = ieo_4stencil[ETA] - j;

            int l0 = mygrid->idx1D(i, j, k);

            for (int jstencil = 0; jstencil < size_of_boundaryStencil; jstencil++) {

              int l1 = leo - mygrid->nXi * jstencil;
              dflux[l0] += multi_fac_4boundaryRight * stencil_boundary[jlocal][jstencil] * func[l1];

            } // jstencil
          } // j
        } // i
    } // mygrid->irank_next[idir_operator][RIGHT]

    break;

  case ZETA:

    // interior scheme
    for (int k = is_4stencil[ZETA]; k <= ie_4stencil[ZETA]; k++)
      for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++)
        for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {

          int l0 = mygrid->idx1D(i, j, k);

          for (int l = -half_size_of_stencil; l <= half_size_of_stencil; l++) {

            int l1 = l0 + mygrid->nXi * mygrid->nEta * l;
            dflux[l0] += stencil[index_abs(l)] * func[l1];

          } // l
        } // i

    // boundary scheme
    // left boundary
    if (mygrid->irank_next[idir_operator][LEFT] == NONE) {
      for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++)
        for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {

          int lso = mygrid->idx1D(i, j, iso_4stencil[ZETA]);
          for (int k = iso_4stencil[ZETA]; k < is_4stencil[ZETA]; k++) {
            int klocal = k - iso_4stencil[ZETA];

            int l0 = mygrid->idx1D(i, j, k);

            for (int kstencil = 0; kstencil < size_of_boundaryStencil; kstencil++) {

              int l1 = lso + mygrid->nXi * mygrid->nEta * kstencil;
              dflux[l0] += stencil_boundary[klocal][kstencil] * func[l1];

            } // kstencil
          } // k
        } // i
    } // mygrid->irank_next[idir_operator][LEFT]

    // right boundary
    if (mygrid->irank_next[idir_operator][RIGHT] == NONE) {
      for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++)
        for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {

          int leo = mygrid->idx1D(i, j, ieo_4stencil[ZETA]);
          for (int k = ieo_4stencil[ZETA]; k > ie_4stencil[ZETA]; k--) {
            int klocal = ieo_4stencil[ZETA] - k;

            int l0 = mygrid->idx1D(i, j, k);

            for (int kstencil = 0; kstencil < size_of_boundaryStencil; kstencil++) {

              int l1 = leo - mygrid->nXi * mygrid->nEta * kstencil;
              dflux[l0] += multi_fac_4boundaryRight * stencil_boundary[klocal][kstencil] * func[l1];

            } // kstencil
          } // k
        } // i
    } // mygrid->irank_next[idir_operator][RIGHT]

    break;

  default:

    mpi::graceful_exit("Unknown direction of taking a derivative.");

    break;

  } // idir_operator

  // correcting derivatives based upon iblank information
  iblank_correction(func, dflux, mygrid, idir_operator);

  return;

} // StandardCentral::apply_spatial_discretization



double StandardCentral::get_multiplicative_factor_4rightBoundary() {

  int multi_fac = DUMMY_DOUBLE;

  switch ( symmetry ) {
  case SYMMETRIC:

    multi_fac = 1.0;
    break;

  case ANTISYMMETRIC:

    multi_fac = -1.0;
    break;

  default:

    mpi::graceful_exit("Unknown configuration of stencil symmetry.");
    break;

  } // symmetry

  return multi_fac;

} // StandardCentral::get_multiplicative_factor_4rightBoundary



void StandardCentral::take_derivative_xi_eta_zeta(double *func, Geometry::StructuredGrid *mygrid, int idir_derivative) {

  // take a derivative of "func" defined on "mygrid" in a direction of "idir_derivative" in a computational coordinate (XI, ETA, ZETA; thus, spacing is unity)
  // the resulting field is stored at dflux[] defined within this namespace spatial

  // initialize
  for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
    dflux[l0] = 0.0;

  // finite difference
  apply_spatial_discretization(func, mygrid, idir_derivative);

  return;

} // StandardCentral::take_derivative_xi_eta_zeta



void StandardCentralFilter::apply_filter(double *func, Geometry::StructuredGrid *mygrid, int idir_filter) {

  // apply a filter on "func" defined on "mygrid" in a direction of "idir_filter" in a computational coordinate (XI, ETA, ZETA; thus, spacing is unity)
  // the resulting field is stored at dflux[] defined within this namespace spatial

  // initialize
  for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
    dflux[l0] = func[l0];

  // filter
  apply_spatial_discretization(func, mygrid, idir_filter);

} // StandardCentralFilter::apply_filter



void StandardCentral::iblank_correction(double *func, double *dfunc, Geometry::StructuredGrid *mygrid, int idir_operator) {

  // my_iblank: iblank value of my cell
  // last_iblank: iblank value of the cell having a one-less index in the direction of interest (idir_operator)
  // e.g. if an i-th cell in the xi direction has a iblank value of "my_iblank", the (i-1)-th cell in the xi direction has "last_iblank"

  int *ijk = new int[DIM_MAX];
  int *iblanks = new int[size_of_boundaryStencil];

  switch ( idir_operator ) {
  case XI:

    for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++) {
      ijk[ZETA] = k;

      for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++) {
        ijk[ETA] = j;

        int l0 = mygrid->idx1D(iso_4stencil[XI], j, k);
        int last_iblank = mygrid->cell[l0].iblank; // initially, last_iblank = my_iblank

        for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {
          ijk[XI] = i;

          l0 = mygrid->idx1D(i, j, k);
          int my_iblank = mygrid->cell[l0].iblank;

          switch ( my_iblank ) {
          case BLANKED: // if I am a hole point (i.e. blanked out, masked)

            if (last_iblank == BLANKED) { // my left cell is a hole as well; blank me out and keep going to the right

              dfunc[l0] = 0.0;

            } // last_iblank
            else if (last_iblank /= BLANKED) { // I am a hole point, but my left cell is not; need to insert right-boundary stencils there

              if (i <= mygrid->is[XI]) { // if my left cell is a ghost, keep going to the right because no need to evaluate derivatives at ghost cells

                // do nothing

              } // i
              else { // insert right-boundary stencils pivoted at my left cell

                // make sure that my grid provides a sufficient number of fluid cells supporting boundary stencils to my left
                check_boundary_stencils_if_all_fluid_points(mygrid, ijk, iblanks, idir_operator, RIGHT);

                // if okay, insert right-boundary stencils
                insert_right_boundary_stencils(mygrid, ijk, func, dfunc, idir_operator);

              } // i
            } // last_iblank

            break;

//          case NOTBLANKED: // I am a regular fluid point on which governing equations are solved
          default:

//            if (last_iblank == NOTBLANKED) { // my left cell is not blanked as well; keep going to the right
//
//              // do nothing
//
//            } // last_iblank
//            else if (last_iblank == BLANKED) {
            if (last_iblank == BLANKED) {

              if (i > mygrid->ie[XI]) { // if I am a ghost, do nothing

                // do nothing

              } // i
              else { // insert left-boundary stencils pivoted at my cell

                // make sure that my grid provides a sufficient number of fluid cells supporting boundary stencils to my right
                check_boundary_stencils_if_all_fluid_points(mygrid, ijk, iblanks, idir_operator, LEFT);

                // if okay, insert left-boundary stencils
                insert_left_boundary_stencils(mygrid, ijk, func, dfunc, idir_operator);

              } // i
            } // last_iblank

            break;

//          default:
//
//            std::cout << ">>> Unknown iblank " << std::setw(4) << my_iblank << " detected at core " << std::setw(4) << mpi::irank << std::endl;
//            mpi::graceful_exit("Unknown iblank.");
//
//            break;

          } // my_iblank

          last_iblank = my_iblank; // update last_iblank and move on to the next cell in the direction idir_operator

        } // i
      } // j
    } // k

    break;

  case ETA:

    for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++) {
      ijk[ZETA] = k;

      for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {
        ijk[XI] = i;

        int l0 = mygrid->idx1D(i, iso_4stencil[ETA], k);
        int last_iblank = mygrid->cell[l0].iblank; // initially, last_iblank = my_iblank

        for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++) {
          ijk[ETA] = j;

          l0 = mygrid->idx1D(i, j, k);
          int my_iblank = mygrid->cell[l0].iblank;

          switch ( my_iblank ) {
          case BLANKED: // if I am a hole point (i.e. blanked out, masked)

            if (last_iblank == BLANKED) { // my left cell is a hole as well; blank me out and keep going to the right

              dfunc[l0] = 0.0;

            } // last_iblank
            else if (last_iblank /= BLANKED) { // I am a hole point, but my left cell is not; need to insert right-boundary stencils there

              if (j <= mygrid->is[ETA]) { // if my left cell is a ghost, keep going to the right because no need to evaluate derivatives at ghost cells

                // do nothing

              } // j
              else { // insert right-boundary stencils pivoted at my left cell

                // make sure that my grid provides a sufficient number of fluid cells supporting boundary stencils to my left
                check_boundary_stencils_if_all_fluid_points(mygrid, ijk, iblanks, idir_operator, RIGHT);

                // if okay, insert right-boundary stencils
                insert_right_boundary_stencils(mygrid, ijk, func, dfunc, idir_operator);

              } // j
            } // last_iblank

            break;

//          case NOTBLANKED: // I am a regular fluid point on which governing equations are solved
          default:

//            if (last_iblank == NOTBLANKED) { // my left cell is not blanked as well; keep going to the right
//
//              // do nothing
//
//            } // last_iblank
//            else if (last_iblank == BLANKED) {
            if (last_iblank == BLANKED) {

              if (j > mygrid->ie[ETA]) { // if I am a ghost, do nothing

                // do nothing

              } // j
              else { // insert left-boundary stencils pivoted at my cell

                // make sure that my grid provides a sufficient number of fluid cells supporting boundary stencils to my right
                check_boundary_stencils_if_all_fluid_points(mygrid, ijk, iblanks, idir_operator, LEFT);

                // if okay, insert left-boundary stencils
                insert_left_boundary_stencils(mygrid, ijk, func, dfunc, idir_operator);

              } // j
            } // last_iblank

            break;

//          default:
//
//            std::cout << ">>> Unknown iblank " << std::setw(4) << my_iblank << " detected at core " << std::setw(4) << mpi::irank << std::endl;
//            mpi::graceful_exit("Unknown iblank.");
//
//            break;

          } // my_iblank

          last_iblank = my_iblank; // update last_iblank and move on to the next cell in the direction idir_operator

        } // j
      } // i
    } // k

    break;

  case ZETA:

    for (int j = iso_4stencil[ETA]; j <= ieo_4stencil[ETA]; j++) {
      ijk[ETA] = j;

      for (int i = iso_4stencil[XI]; i <= ieo_4stencil[XI]; i++) {
        ijk[XI] = i;

        int l0 = mygrid->idx1D(i, j, iso_4stencil[ZETA]);
        int last_iblank = mygrid->cell[l0].iblank; // initially, last_iblank = my_iblank

        for (int k = iso_4stencil[ZETA]; k <= ieo_4stencil[ZETA]; k++) {
          ijk[ZETA] = k;

          l0 = mygrid->idx1D(i, j, k);
          int my_iblank = mygrid->cell[l0].iblank;

          switch ( my_iblank ) {
          case BLANKED: // if I am a hole point (i.e. blanked out, masked)

            if (last_iblank == BLANKED) { // my left cell is a hole as well; blank me out and keep going to the right

              dfunc[l0] = 0.0;

            } // last_iblank
            else if (last_iblank /= BLANKED) { // I am a hole point, but my left cell is not; need to insert right-boundary stencils there

              if (k <= mygrid->is[ZETA]) { // if my left cell is a ghost, keep going to the right because no need to evaluate derivatives at ghost cells

                // do nothing

              } // k
              else { // insert right-boundary stencils pivoted at my left cell

                // make sure that my grid provides a sufficient number of fluid cells supporting boundary stencils to my left
                check_boundary_stencils_if_all_fluid_points(mygrid, ijk, iblanks, idir_operator, RIGHT);

                // if okay, insert right-boundary stencils
                insert_right_boundary_stencils(mygrid, ijk, func, dfunc, idir_operator);

              } // k
            } // last_iblank

            break;

//          case NOTBLANKED: // I am a regular fluid point on which governing equations are solved
          default:

//            if (last_iblank == NOTBLANKED) { // my left cell is not blanked as well; keep going to the right
//
//              // do nothing
//
//            } // last_iblank
//            else if (last_iblank == BLANKED) {
            if (last_iblank == BLANKED) {

              if (k > mygrid->ie[ZETA]) { // if I am a ghost, do nothing

                // do nothing

              } // k
              else { // insert left-boundary stencils pivoted at my cell

                // make sure that my grid provides a sufficient number of fluid cells supporting boundary stencils to my right
                check_boundary_stencils_if_all_fluid_points(mygrid, ijk, iblanks, idir_operator, LEFT);

                // if okay, insert left-boundary stencils
                insert_left_boundary_stencils(mygrid, ijk, func, dfunc, idir_operator);

              } // k
            } // last_iblank

            break;

//          default:
//
//            std::cout << ">>> Unknown iblank " << std::setw(4) << my_iblank << " detected at core " << std::setw(4) << mpi::irank << std::endl;
//            mpi::graceful_exit("Unknown iblank.");
//
//            break;

          } // my_iblank

          last_iblank = my_iblank; // update last_iblank and move on to the next cell in the direction idir_operator

        } // k
      } // j
    } // i

    break;

  default:

    mpi::graceful_exit("Unknown direction of taking a derivative.");

    break;

  } // idir_operator

  DEALLOCATE_1DPTR(ijk);
  DEALLOCATE_1DPTR(iblanks);

  return;

} // StandardCentral::iblank_correction



void StandardCentral::insert_right_boundary_stencils(Geometry::StructuredGrid *mygrid, int *ijk_in_grid, double *func, double *dfunc, int idir_operator) {

  // ijk_in_grid is assumed a hole point having a fluid point to its left;
  // thus, at either (i-1,j,k), (i,j-1,k), or (i,j,k-1), boundary stencils are inserted to the left

  int shift_in_stencil[DIM_MAX];
  for (int idir = XI; idir < DIM_MAX; idir++)
    shift_in_stencil[idir] = 0;

  int ijk_work[DIM_MAX];

  int shift_in_cells[DIM_MAX];
  shift_in_cells[XI] = 1; // -1 because boundary stencils start growing (to the left) at either i-1 or j-1 or k-1
  shift_in_cells[ETA] = mygrid->num_ocells_dir[XI];
  shift_in_cells[ZETA] = mygrid->num_ocells_dir[XI] * mygrid->num_ocells_dir[ETA];

  // get the reference 1-D index for the boundary (where the stencil is most biased to the left)
  shift_in_stencil[idir_operator] = -1;
  for (int idir = XI; idir < DIM_MAX; idir++)
    ijk_work[idir] = ijk_in_grid[idir] + shift_in_stencil[idir];
  int lref = mygrid->idx1D(ijk_work[XI], ijk_work[ETA], ijk_work[ZETA]);

  for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) {

    shift_in_stencil[idir_operator] = -1 - iboundary; // -1 because boundary stencils start growing (to the left) at either i-1 or j-1 or k-1
    for (int idir = XI; idir < DIM_MAX; idir++)
      ijk_work[idir] = ijk_in_grid[idir] + shift_in_stencil[idir]; // grid-level ijk-index of this boundary cell at iboundary from the right

    if ((ijk_work[idir_operator] - mygrid->is[idir_operator]) * (ijk_work[idir_operator] - mygrid->ie[idir_operator]) <= 0) { // only non-ghost cells are updated

      int l1 = mygrid->idx1D(ijk_work[XI], ijk_work[ETA], ijk_work[ZETA]); // grid-level 1-D index of this boundary cell at iboundary from the right

      dfunc[l1] = 0.0;
      for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++) {

        int l2 = lref - istencil * shift_in_cells[idir_operator];
        dfunc[l1] += multi_fac_4boundaryRight * stencil_boundary[iboundary][istencil] * func[l2];

      } // istencil
    } // (ijk_work[idir_operator] - mygrid->is[idir_operator])
  } // iboundary

  return;

} // StandardCentral::insert_right_boundary_stencils



void StandardCentral::insert_left_boundary_stencils(Geometry::StructuredGrid *mygrid, int *ijk_in_grid, double *func, double *dfunc, int idir_operator) {

  // ijk_in_grid is assumed a fluid point having a hole point to its left;
  // thus, at (i,j,k), boundary stencils are inserted to the right

  int shift_in_stencil[DIM_MAX];
  for (int idir = XI; idir < DIM_MAX; idir++)
    shift_in_stencil[idir] = 0;

  int ijk_work[DIM_MAX];

  int shift_in_cells[DIM_MAX];
  shift_in_cells[XI] = 1;
  shift_in_cells[ETA] = mygrid->num_ocells_dir[XI];
  shift_in_cells[ZETA] = mygrid->num_ocells_dir[XI] * mygrid->num_ocells_dir[ETA];

  // get the reference 1-D index for the boundary (where the stencil is most biased to the right)
  int lref = mygrid->idx1D(ijk_in_grid[XI], ijk_in_grid[ETA], ijk_in_grid[ZETA]);

  for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) {

    shift_in_stencil[idir_operator] = iboundary;
    for (int idir = XI; idir < DIM_MAX; idir++)
      ijk_work[idir] = ijk_in_grid[idir] + shift_in_stencil[idir]; // grid-level ijk-index of this boundary cell at iboundary from the left

    if ((ijk_work[idir_operator] - mygrid->is[idir_operator]) * (ijk_work[idir_operator] - mygrid->ie[idir_operator]) <= 0) { // only non-ghost cells are updated

      int l1 = mygrid->idx1D(ijk_work[XI], ijk_work[ETA], ijk_work[ZETA]); // grid-level 1-D index of this boundary cell at iboundary from the left

      dfunc[l1] = 0.0;
      for (int istencil = 0; istencil < size_of_boundaryStencil; istencil++) {

        int l2 = lref + istencil * shift_in_cells[idir_operator];
        dfunc[l1] += stencil_boundary[iboundary][istencil] * func[l2];

      } // istencil
    } // (ijk_work[idir_operator] - mygrid->is[idir_operator])
  } // iboundary

  return;

} // StandardCentral::insert_left_boundary_stencils



void StandardCentral::check_boundary_stencils_if_all_fluid_points(Geometry::StructuredGrid *mygrid, int *ijk_in_grid, int *iblanks, int idir_operator, int which_end) {

  int shift_in_stencil[DIM_MAX];
  for (int idir = XI; idir < DIM_MAX; idir++)
    shift_in_stencil[idir] = 0;

  int ijk_work[DIM_MAX];

  switch ( which_end ) {
  case LEFT: // it is assumed that ijk_in_grid or (i,j,k) is a first fluid point after a number of hole points
             // thus, ijk_in_grid or (i,j,k) is where boundary stencils start growing (to the right)

    for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) {

      shift_in_stencil[idir_operator] = iboundary;
      for (int idir = XI; idir < DIM_MAX; idir++)
        ijk_work[idir] = ijk_in_grid[idir] + shift_in_stencil[idir]; // grid-level ijk-index of this boundary cell at iboundary from the left

      int l1 = mygrid->idx1D(ijk_work[XI], ijk_work[ETA], ijk_work[ZETA]); // grid-level 1-D index of this boundary cell at iboundary from the left
      iblanks[iboundary] = std::max( mygrid->cell[l1].iblank, 0 ); // max. operation so that interpolated cells having a negative value of iblank can also serve as stencil members

    } // iboundary

    break;

  case RIGHT: // it is assumed that ijk_in_grid or (i,j,k) is a first hole point after a number of fluid points
              // thus, either (i-1,j,k) or (i,j-1,k) or (i,j,k-1) is where boundary stencils start growing (to the left)

    for (int iboundary = FIRST; iboundary < num_of_boundaryCells; iboundary++) {

      shift_in_stencil[idir_operator] = -1 - iboundary; // -1 because boundary stencils start growing (to the left) at either i-1 or j-1 or k-1
      for (int idir = XI; idir < DIM_MAX; idir++)
        ijk_work[idir] = ijk_in_grid[idir] + shift_in_stencil[idir]; // grid-level ijk-index of this boundary cell at iboundary from the right

      int l1 = mygrid->idx1D(ijk_work[XI], ijk_work[ETA], ijk_work[ZETA]); // grid-level 1-D index of this boundary cell at iboundary from the right
      iblanks[iboundary] = std::max( mygrid->cell[l1].iblank, 0 ); // max. operation so that interpolated cells having a negative value of iblank can also serve as stencil members

    } // iboundary

    break;

  } // which_end

  // ensure such that there are no hole points within the inserted boundary stencils; 
  // otherwise, there are too few fluid points between two adjacent holes
  assert ( math_algebra::sum_array_integer(iblanks, num_of_boundaryCells) == NOTBLANKED );

  return;

} // StandardCentral::check_boundary_stencils_if_all_fluid_points



void StandardCentral::take_derivative_xyz(double *func, Geometry::StructuredGrid *mygrid, int idir_derivative) {

  // takes a derivative of "func" defined on "mygrid" in a direction of "idir_derivative" in a Cartesian coordinate (x, y, z)

  double *dflux_sum = new double[mygrid->num_ocells];
  for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
    dflux_sum[l0] = 0.0;

  for (int l = XI; l < num_dim; l++) {

    this->take_derivative_xi_eta_zeta(func, mygrid, l); // d/dxi_l

    for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
      dflux_sum[l0] += mygrid->cell[l0].metrics[l][idir_derivative] * dflux[l0]; // \hat{dxi_l / dx_i} * d/dxi_l

  } // l

  for (int l0 = 0; l0 < mygrid->num_ocells; l0++)
    dflux[l0] = dflux_sum[l0] * mygrid->cell[l0].Jac;

  DEALLOCATE_1DPTR(dflux_sum);

  return;

} // StandardCentral::take_derivative_xyz



void initialize(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate) {

  std::string name_scheme;

  num_dim = mystate->num_dim;
  num_vars = mystate->num_vars_sol;
  num_samples = mystate->num_samples;

  // finite difference
  spatial_scheme = myinput->spatial_scheme;
  if (spatial_scheme == "STANDARD_CENTRAL") {

    cds.initialize(myinput, mygrid);
    name_scheme = cds.name_scheme;

  } // spatial_scheme
  else if (spatial_scheme == "FDO11P") {

    cdo.initialize(myinput, mygrid);
    name_scheme = cdo.name_scheme;

  } // spatial_scheme
  else
    mpi::graceful_exit("Unknown spatial scheme.");
  //
  mpi::wait_allothers("Spatial discretization is " + name_scheme + " and is initialized.");

//  std::cout << "Rank: " << std::setw(4) << mpi::irank << ", cds operator is-ie in XI: " << cds.is_4stencil[XI]   << "-"  << cds.ie_4stencil[XI]
//                                                      << ", in ETA: "                   << cds.is_4stencil[ETA]  << "-"  << cds.ie_4stencil[ETA]
//                                                      << ", in ZETA: "                  << cds.is_4stencil[ZETA] << "-"  << cds.ie_4stencil[ZETA] << std::endl;

  // filter
  do_filter = myinput->do_filter;
  filter_scheme = myinput->filter_scheme;
  if (do_filter == TRUE) {
    if (filter_scheme == "STANDARD_CENTRAL") {

      cfs.initialize(myinput, mygrid);
      name_scheme = cfs.name_scheme;

    } // filter_scheme
    else if (filter_scheme == "SFO11P") {

      cfo.initialize(myinput, mygrid);
      name_scheme = cfo.name_scheme;

    } // filter_scheme
    else
      mpi::graceful_exit("Unknown filter scheme.");
    //
    mpi::wait_allothers("Filter is " + name_scheme + " and is initialized.");

//    std::cout << "Rank: " << std::setw(4) << mpi::irank << ", cfs operator is-ie in XI: " << cfs.is_4stencil[XI]   << "-"  << cfs.ie_4stencil[XI]
//                                                        << ", in ETA: "                   << cfs.is_4stencil[ETA]  << "-"  << cfs.ie_4stencil[ETA]
//                                                        << ", in ZETA: "                  << cfs.is_4stencil[ZETA] << "-"  << cfs.ie_4stencil[ZETA] << std::endl;

  } // do_filter

  // boundary condition
  if (spatial_scheme == "STANDARD_CENTRAL")
    bc::initialize_bc(num_vars, cds.num_of_boundaryCells, cds.size_of_boundaryStencil, cds.stencil_boundary);

  else if (spatial_scheme == "FDO11P")
    bc::initialize_bc(num_vars, cdo.num_of_boundaryCells, cdo.size_of_boundaryStencil, cdo.stencil_boundary);

  else
    mpi::graceful_exit("Unknown spatial scheme.");

  // work storage for spatial discretization
  flux = new double[num_samples];
  dflux = new double[num_samples];

  return;

} // initialize



void finalize() {

  DEALLOCATE_1DPTR(flux);
  DEALLOCATE_1DPTR(dflux);

  return;

} // finalize



void take_derivative_xi_eta_zeta(double *func, Geometry::StructuredGrid *mygrid, int idir_derivative) {

  if (spatial_scheme == "STANDARD_CENTRAL")
    cds.take_derivative_xi_eta_zeta(func, mygrid, idir_derivative);

  else if (spatial_scheme == "FDO11P")
    cdo.take_derivative_xi_eta_zeta(func, mygrid, idir_derivative);

  else
    mpi::graceful_exit("Unknown spatial scheme.");

  return;

} // take_derivative_xi_eta_zeta



void take_derivative_xyz(double *func, Geometry::StructuredGrid *mygrid, int idir_derivative) {

  if (spatial_scheme == "STANDARD_CENTRAL")
    cds.take_derivative_xyz(func, mygrid, idir_derivative);

  else if (spatial_scheme == "FDO11P")
    cdo.take_derivative_xyz(func, mygrid, idir_derivative);

  else
    mpi::graceful_exit("Unknown spatial scheme.");

  return;

} // take_derivative_xyz



void zero_out_flux() {

  for (int l0 = 0; l0 < num_samples; l0++)
    flux[l0] = 0.0;

  return;

}


void compute_RHS(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs, double time) {

  // zero out RHS
  for (int ivar = 0; ivar < num_vars; ivar++)
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[ivar])[l0] = 0.0;

  // physical model
  if (myinput->model_pde == "LINEAR_ACOUSTICS")
    compute_RHS_acoustics(myinput, mygrid, mystate, y, rhs);

  else if (myinput->model_pde == "LEE")
    compute_RHS_linearizedEuler(myinput, mygrid, mystate, y, rhs);

  else if (myinput->model_pde == "LEE_SCALAR") {
    compute_RHS_linearizedEuler(myinput, mygrid, mystate, y, rhs);
    compute_RHS_linearizedEuler_scalar(myinput, mygrid, mystate, y, rhs);

  } // myinput->model_pde
  else if (myinput->model_pde == "LEE_MIXFRAC_CONSTGAMMA") {
    compute_RHS_linearizedEuler(myinput, mygrid, mystate, y, rhs);
    compute_RHS_linearizedEuler_scalar(myinput, mygrid, mystate, y, rhs);

  } // myinput->model_pde
  else
    mpi::graceful_exit("PHYSICAL_MODEL = " + myinput->model_pde + " is unknown, so the RHS of its governing equations cannot be evaluated.");

  // compute the RHS source terms
  source::add_RHSSourceTerms(myinput, mygrid, mystate, y, rhs, time);

  return;

} // compute_RHS



void compute_RHS_acoustics(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs) {

  double rho_0 = 1.0; // ambient density
  double p_0 = 1.0 / mystate->gamma_specificheat; // ratio of specific heats, gamma = C_p / C_v
  double c_0 = 1.0; // ambient speed of sound

  // velocity fluctuations
  for (int idir_vel = XDIR; idir_vel < num_dim; idir_vel++) { // idir_vel denotes an index for a velocity equation
    for (int idir_drv = XI; idir_drv < num_dim; idir_drv++) { // idir_drv denotes a direction (xi, eta, zeta) in which you are taking a derivative

      for (int l0 = 0; l0 < num_samples; l0++)
        flux[l0] = mygrid->cell[l0].metrics[idir_drv][idir_vel] * (y[IVAR_P])[l0];

      take_derivative_xi_eta_zeta(flux, mygrid, idir_drv);

      for (int l0 = 0; l0 < num_samples; l0++)
        (rhs[IVAR_UX + idir_vel])[l0] -= dflux[l0];

    } // idir_drv
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_UX + idir_vel])[l0] *= mygrid->cell[l0].Jac / rho_0;

  } // idir_vel

  // pressure fluctuation
  for (int idir_drv = XI; idir_drv < num_dim; idir_drv++) {

    zero_out_flux();
    for (int idir_vel = XDIR; idir_vel < num_dim; idir_vel++)
      for (int l0 = 0; l0 < num_samples; l0++)
        flux[l0] += mygrid->cell[l0].metrics[idir_drv][idir_vel] * (y[IVAR_UX + idir_vel])[l0];

    take_derivative_xi_eta_zeta(flux, mygrid, idir_drv);

    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_P])[l0] -= dflux[l0];

  } // idir_drv
  for (int l0 = 0; l0 < num_samples; l0++)
    (rhs[IVAR_P])[l0] *= mystate->gamma_specificheat * p_0 * mygrid->cell[l0].Jac;

  // additional terms due to axisymmetry
  if (myinput->axisym == TRUE) {

    // in the following, Jacobian is not needed since it is already multiplied with derivatives above
    // (the additional term arising from axisymmetry does not contain any derivative)
    // also, 1/r may cause a numerical issue but since this is axisymmetric, a boundary condition for pressure will be imposed at r = 0
    // thus, it is ok just to add a small number to the denominator to avoid a singularity for now
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_P])[l0] -= mystate->gamma_specificheat * p_0 * (y[IVAR_UR])[l0] / (mygrid->cell[l0].xyz[RDIR] + DUMMY_SMALL);

  } // myinput->axisym

  return;

} // compute_RHS_acoustics



void compute_RHS_linearizedEuler(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs) {

  // this RHS-computing routine calculates xyz derivatives and does not assume a strong conservation form

  // entropy fluctuation
  for (int idir_drv = XDIR; idir_drv < num_dim; idir_drv++) {

    // -\bar{u}_i \frac{\partial s^\prime}{\partial x_i}
    take_derivative_xyz(y[IVAR_S], mygrid, idir_drv);
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_S])[l0] -= (mystate->sol_mean[IVAR_UX + idir_drv])[l0] * dflux[l0];

    // -u^\prime_i \frac{\partial \bar{s}}{\partial x_i}
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_S])[l0] -= (y[IVAR_UX + idir_drv])[l0] * (mystate->sol_meanGradient[mystate->ivar1D(IVAR_S, idir_drv)])[l0];

  } // idir_drv

  // velocity fluctuations
  for (int idir_vel = XDIR; idir_vel < num_dim; idir_vel++) { // idir_vel denotes an index for a velocity equation
    for (int idir_drv = XDIR; idir_drv < num_dim; idir_drv++) { // idir_drv denotes a direction (x, y, z) in which you are taking a derivative

      // -\bar{u}_j \frac{\partial u^\prime_i}{\partial x_j}
      take_derivative_xyz(y[IVAR_UX + idir_vel], mygrid, idir_drv);
      for (int l0 = 0; l0 < num_samples; l0++)
        (rhs[IVAR_UX + idir_vel])[l0] -= (mystate->sol_mean[IVAR_UX + idir_drv])[l0] * dflux[l0];

      // -(u^\prime_j + \frac{\rho^\prime}{\bar{\rho}} \bar{u}_j) \frac{\partial \bar{u}_i}{\partial x_j}
      for (int l0 = 0; l0 < num_samples; l0++)
        (rhs[IVAR_UX + idir_vel])[l0] -= ((y[IVAR_UX + idir_drv])[l0]
                                       + (mystate->sol_aux[IAUX_RHO])[l0] / (mystate->sol_aux[IAUX_RHO_MEAN])[l0] * (mystate->sol_mean[IVAR_UX + idir_drv])[l0])
                                       * (mystate->sol_meanGradient[mystate->ivar1D(IVAR_UX + idir_vel, idir_drv)])[l0];

    } // idir_drv

    // -\frac{1}{\bar{\rho}} \frac{\partial p^\prime}{\partial x_i}
    take_derivative_xyz(y[IVAR_P], mygrid, idir_vel);
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_UX + idir_vel])[l0] -= dflux[l0] / (mystate->sol_aux[IAUX_RHO_MEAN])[l0];

  } // idir_vel

  // pressure fluctuation
  for (int idir_drv = XDIR; idir_drv < num_dim; idir_drv++) {

    // -\bar{u}_i \frac{\partial p^\prime}{\partial x_i}
    take_derivative_xyz(y[IVAR_P], mygrid, idir_drv);
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_P])[l0] -= (mystate->sol_mean[IVAR_UX + idir_drv])[l0] * dflux[l0];

    // -u^\prime_i \frac{\partial \bar{p}}{\partial x_i}
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_P])[l0] -= (y[IVAR_UX + idir_drv])[l0] * (mystate->sol_meanGradient[mystate->ivar1D(IVAR_P, idir_drv)])[l0];

    // -\gamma \bar{p} \frac{\partial u^\prime_i}{\partial x_i}
    take_derivative_xyz(y[IVAR_UX + idir_drv], mygrid, idir_drv);
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_P])[l0] -= mystate->gamma_specificheat * (mystate->sol_mean[IVAR_P])[l0] * dflux[l0];

    // -\gamma p^\prime \frac{\partial \bar{u}_i}{\partial x_i}
    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_P])[l0] -= mystate->gamma_specificheat * (y[IVAR_P])[l0] * (mystate->sol_meanGradient[mystate->ivar1D(IVAR_UX + idir_drv, idir_drv)])[l0];

  } // idir_drv

  // additional terms due to axisymmetry
  if (myinput->axisym == TRUE) {

    for (int l0 = 0; l0 < num_samples; l0++)
      (rhs[IVAR_P])[l0] -= mystate->gamma_specificheat
                        * ((mystate->sol_mean[IVAR_P])[l0] * (y[IVAR_UR])[l0] + (y[IVAR_P])[l0] * (mystate->sol_mean[IVAR_UR])[l0])
                        / (mygrid->cell[l0].xyz[RDIR] + DUMMY_SMALL);

  } // myinput->axisym

  return;

} // compute_RHS_linearizedEuler



void compute_RHS_linearizedEuler_scalar(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, double **y, double **rhs) {

  // this RHS-computing routine calculates xyz derivatives and does not assume a strong conservation form

  int ivar_shift = IVAR_P + 1;

  // scalar fluctuation
  for (int ivar = 0; ivar < myinput->num_scalar; ivar++) {
    for (int idir_drv = XDIR; idir_drv < num_dim; idir_drv++) {

      // -\bar{u}_i \frac{\partial Z^\prime}{\partial x_i}
      take_derivative_xyz(y[ivar_shift+ivar], mygrid, idir_drv);
      for (int l0 = 0; l0 < num_samples; l0++)
        (rhs[ivar_shift+ivar])[l0] -= (mystate->sol_mean[IVAR_UX + idir_drv])[l0] * dflux[l0];

      // -u^\prime_i \frac{\partial \bar{Z}}{\partial x_i}
      for (int l0 = 0; l0 < num_samples; l0++)
        (rhs[ivar_shift+ivar])[l0] -= (y[IVAR_UX + idir_drv])[l0] * (mystate->sol_meanGradient[mystate->ivar1D(ivar_shift+ivar, idir_drv)])[l0];

    } // idir_drv
  } // ivar

  return;

} // compute_RHS_linearizedEuler_scalar



void update_boundary(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate, int num_vars_in, double **y, double time) {

  // physical boundary conditions
  bc::apply_BC(myinput, mygrid, mystate, num_vars_in, y, time);

  // overset grid interpolation
  if (myinput->do_overset == TRUE) {

    // since solutions are time advanced BEFORE applying boundary conditions and ghost cells are not time advanced, 
    // we need to additionally update ghost-cell data in every direction
    for (int idir = XI; idir < num_dim; idir++)
      update_ghostcell_data(y, num_vars_in, mygrid, idir);

    // the following implementation does message passing num_vars_in times more
    // the function itself is 30% slower than above
//    for (int ivar = 0; ivar < num_vars_in; ivar++)
//      for (int idir = XI; idir < num_dim; idir++)
//        update_ghostcell_data(y[ivar], mygrid, idir);

    overset::interpolate(mygrid, y);

  } // myinput->do_overset

  return;

} // update_boundary



void apply_filter(int num_vars_in, int num_samples_in, double **y, Geometry::StructuredGrid *mygrid) {

  int dir_filter[DIM_MAX]; // direction in which filter is applied
  dir_filter[FIRST] = XI;
  dir_filter[SECOND] = ETA;
  dir_filter[THIRD] = ZETA;
  double filter_blend, filter_blend_remainder;

  double *func0 = new double[num_samples_in];
  double *func = new double[num_samples_in];

  // filter each variable in each direction sequentially
  for (int ivar = 0; ivar < num_vars_in; ivar++) {

    for (int l0 = 0; l0 < num_samples_in; l0++)
      func0[l0] = (y[ivar])[l0]; // keep the unfilterd values to blend with the filtered values
    for (int l0 = 0; l0 < num_samples_in; l0++)
      func[l0] = func0[l0];

    for (int idir = FIRST; idir < num_dim; idir++) {

      int idir_filter = dir_filter[idir];

      if (filter_scheme == "STANDARD_CENTRAL") {

        filter_blend = cfs.filter_blend;
        cfs.apply_filter(func, mygrid, idir_filter);

      } // filter_scheme
      else if (filter_scheme == "SFO11P") {

        filter_blend = cfo.filter_blend;
        cfo.apply_filter(func, mygrid, idir_filter);

      } // filter_scheme
      else
        mpi::graceful_exit("Unknown filter scheme.");

      for (int l0 = 0; l0 < num_samples_in; l0++)
        func[l0] = dflux[l0];

    } // idir

    // blend unfiltered and filterd values
    filter_blend_remainder = 1.0 - filter_blend;
    for (int l0 = 0; l0 < num_samples_in; l0++)
      (y[ivar])[l0] = filter_blend_remainder * func0[l0] + filter_blend * func[l0];

  } // ivar

  DEALLOCATE_1DPTR(func0);
  DEALLOCATE_1DPTR(func);

  return;

} // apply_filter



void precompute_something(UserInput *myinput, Geometry::StructuredGrid *mygrid, State *mystate) {

  // xyz gradients of base (or mean) state
  if (mystate->model_pde == "LEE" || 
      mystate->model_pde == "LEE_SCALAR" ||
      mystate->model_pde == "LEE_MIXFRAC_CONSTGAMMA") {

    for (int ivar = 0; ivar < mystate->num_vars_mean; ivar++) {
      for (int idir_drv = XDIR; idir_drv < num_dim; idir_drv++) {

        take_derivative_xyz(mystate->sol_mean[ivar], mygrid, idir_drv);
        for (int l0 = 0; l0 < mystate->num_samples; l0++)
          (mystate->sol_meanGradient[mystate->ivar1D(ivar, idir_drv)])[l0] = dflux[l0];

      } // idir_drv
    } // ivar

  } // mystate->model_pde

  return;

} // precompute_something

} // spatial
