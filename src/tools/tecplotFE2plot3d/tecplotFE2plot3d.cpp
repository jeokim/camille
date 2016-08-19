#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "../../core/simulation.h"
#include "../../io/tecplot.h"
#include "../../math/matrix.h"

const double eps_smallAngle = 0.05; //0.0;
const double eps_proximity = 0.0; //pow(10.0, -6);

int main(int argc, char * argv[]) {

  // global variables and arrays
  int num_zones_source;
  int num_dim_source;
  int num_vars_source;
  std::string vars_source;

  double xyz_min[DIM_MAX], xyz_max[DIM_MAX];

  double xyz_receiver[DIM_MAX];
  int *nodeID_donor;
  double **xyz_donor;
  double **var_donor;

  double vec0[DIM_MAX], vec1[DIM_MAX];

  // simplicity
  UserInput *myinput = &(simulation::myinput);
  Geometry::StructuredBlock *&block = simulation::block;
  Geometry::StructuredGrid *&mygrid = simulation::grid; // "Geometry::StructuredGrid *mygrid = simulation::grid;" doesn't work 
                                                        // since mygrid gets only a copy of "simulation::grid" although the syntax is 
                                                        // consistent (a pointer copied to a pointer); to access the original 
                                                        // "simulation::grid", it is necessary to first make the receiver "mygrid" as 
                                                        // a reference (by putting & in front of "mygrid"); then, to be consistent (a 
                                                        // pointer on LHS and a pointer on RHS), place * in its front
  State *mystate = &(simulation::state);



  // since the interpolation is done using multiple processors and filtering is applied, the
  // entire initialization routine is called
  simulation::initialize(argc, argv);
  if (myinput->num_filters_interpolation > 0 && myinput->do_overset == TRUE)
    mpi::graceful_exit("Somehow, filtering interpolated solution causes problems at overset boundaries; turn off the overset only for interpolation.");



  // some temporary hacks
  if (mpi::irank == 0)
    std::cout << std::endl;
  MESSAGE_STDOUT("Tecplot--PLOT3D interpolation is currently limited to 2D-2D cases (axisymmetric case included).");
  MESSAGE_STDOUT("Tecplot ASCII file should contain absolute pressure, density, axial (or streamwise) velocity, radial (or transverse) velocity.");
  MESSAGE_STDOUT("Alternatively, Tecplot ASCII file can contain absolute pressure, axial (or streamwise) velocity, radial (or transverse) velocity, specific entropy.");
  MESSAGE_STDOUT("Tecplot ASCII file written in FE POINT data-packing is supported only; so write it consistently.");
  if (myinput->num_dim != 2)
    mpi::graceful_exit("Only 2D-2D interpolation is supported.");



  // initialize the interpolation-source data
  num_zones_source = myinput->num_zones_interpSource;
  num_dim_source = myinput->num_dim_interpSource;
  num_vars_source = myinput->num_vars_interpSource;
  vars_source = myinput->vars_interpSource;
  t_TecplotFEZone *TecplotFEZone = new t_TecplotFEZone[num_zones_source];
  for (int izone = 0; izone < num_zones_source; izone++) {

    TecplotFEZone[izone].zoneID = izone;
    TecplotFEZone[izone].num_dim = num_dim_source;
    TecplotFEZone[izone].num_vars = num_vars_source;
    TecplotFEZone[izone].num_vertices_per_cv = 0; // set below
    TecplotFEZone[izone].num_nodes = 0; // set below
    TecplotFEZone[izone].num_cvs = 0; // set below

  } // izone



  // read a Tecplot ASCII file
  // for formatted file-reading, use fgets and fscanf
  FILE *pFile = fopen(cstr_to_constchar(myinput->file_tecplot_ASCII), "r");
  if (pFile == NULL)
    mpi::graceful_exit("The Tecplot ASCII file to interpolate does not exist.");
  //
  tecplot::read_ASCII_header_FE(pFile, num_dim_source, num_vars_source);
  for (int izone = 0; izone < num_zones_source; izone++) {

    // simplicity
    int num_nodes_Tecplot = 0;
    int num_cvs_Tecplot = 0;

    if (mpi::irank == 0)
      std::cout << ">>> Zone " << std::setw(3) << izone << ": " << std::endl;

    TecplotFEZone[izone].num_vertices_per_cv = tecplot::read_ASCII_zoneHeader_FE(pFile, num_dim_source, num_vars_source, num_nodes_Tecplot, num_cvs_Tecplot);
    TecplotFEZone[izone].num_nodes = num_nodes_Tecplot;
    TecplotFEZone[izone].num_cvs = num_cvs_Tecplot;

    tecplot::read_ASCII_POINT(pFile, num_dim_source, num_vars_source, num_nodes_Tecplot, num_cvs_Tecplot, 
                              TecplotFEZone[izone].xyz, 
                              TecplotFEZone[izone].var, 
                              TecplotFEZone[izone].connectivity, 
                              TecplotFEZone[izone].num_vertices_per_cv);

  } // izone
  fclose(pFile);
  mpi::wait_allothers("An interpolation-source file is read and stored.");



  // re-scale the xyz locations of the source data
  double inv_scale_xyz = 1.0 / myinput->scale_xyz;
  for (int izone = 0; izone < num_zones_source; izone++)
    for (int inode = 0; inode < TecplotFEZone[izone].num_nodes; inode++)
      for (int idim = XDIR; idim < num_dim_source; idim++)
        (TecplotFEZone[izone].xyz[inode])[idim] *= inv_scale_xyz;
  // find out and store the minimum and maximum locations
  for (int idim = XDIR; idim < DIM_MAX; idim++) {
    xyz_min[idim] = DUMMY_LARGE;
    xyz_max[idim] = -DUMMY_LARGE;
  } // idim
  for (int izone = 0; izone < num_zones_source; izone++)
    for (int inode = 0; inode < TecplotFEZone[izone].num_nodes; inode++)
      for (int idim = XDIR; idim < num_dim_source; idim++) {
        xyz_min[idim] = std::min(xyz_min[idim], (TecplotFEZone[izone].xyz[inode])[idim]);
        xyz_max[idim] = std::max(xyz_max[idim], (TecplotFEZone[izone].xyz[inode])[idim]);
      } // idim
//  if (mpi::irank == 0) {
//    std::cout << xyz_min[XDIR] << " to " << xyz_max[XDIR] << std::endl;
//    std::cout << xyz_min[YDIR] << " to " << xyz_max[YDIR] << std::endl;
//  } // mpi::irank
  mpi::wait_allothers("The spatial locations of the interpolation-source data are re-scaled.");



  // post-process the solution to be interpolated
  double gamma;
  for (int izone = FIRST; izone < num_zones_source; izone++) {

    // simplicity
    int num_nodes_Tecplot = TecplotFEZone[izone].num_nodes;
    double **var_source = TecplotFEZone[izone].var;

    if (myinput->model_pde == "LINEAR_EULER") {

      // re-scale the solution to be interpolated
      if (vars_source == "PRHOUXUR" || vars_source == "PRHOUXUY") {
        if (izone == FIRST)
          MESSAGE_STDOUT("P, RHO, U_i are assumed to be in the interpolation-source file in that order.");
        double inv_scale_p = 1.0 / myinput->scale_p;
        double inv_scale_rho = 1.0 / myinput->scale_rho;
        double inv_scale_u = 1.0 / myinput->scale_u;
        for (int inode = 0; inode < num_nodes_Tecplot; inode++) {

          var_source[inode][FIRST] *= inv_scale_p; // P
          var_source[inode][SECOND] *= inv_scale_rho; // RHO
          var_source[inode][THIRD] *= inv_scale_u; // U_X
          var_source[inode][FOURTH] *= inv_scale_u; // U_Y

        } // inode
      } // vars_source
      else if (vars_source == "PUXURS" || vars_source == "PUXUYS") {
        if (izone == FIRST)
          MESSAGE_STDOUT("P, U_i, specific entropy are assumed to be in the interpolation-source file in that order.");
        double inv_scale_p = 1.0 / myinput->scale_p;
        double inv_scale_u = 1.0 / myinput->scale_u;
        double inv_scale_s = 1.0 / myinput->scale_s;
        for (int inode = 0; inode < num_nodes_Tecplot; inode++) {

          var_source[inode][FIRST] *= inv_scale_p; // P
          var_source[inode][SECOND] *= inv_scale_u; // U_X
          var_source[inode][THIRD] *= inv_scale_u; // U_Y
          var_source[inode][FOURTH] *= inv_scale_s; // s (specific entropy)

        } // inode
      } // vars_source
      if (izone == num_zones_source - 1)
        mpi::wait_allothers("The solution to be interpolated is non-dimensionalized.");

      gamma = myinput->gamma_specificheat;
      if (vars_source == "PRHOUXUR" || vars_source == "PRHOUXUY") {
        for (int inode = 0; inode < num_nodes_Tecplot; inode++) {

          // convert {P, RHO, U_i} into {S, U_i, P}
          double entropy = (log(gamma * var_source[inode][FIRST]) - gamma * log(var_source[inode][SECOND])) / gamma;
          double pressure = var_source[inode][FIRST];
          double ux = var_source[inode][THIRD];
          double uy = var_source[inode][FOURTH];

          var_source[inode][FIRST] = entropy;
          var_source[inode][SECOND] = ux;
          var_source[inode][THIRD] = uy;
          var_source[inode][FOURTH] = pressure;

        } // inode
      } // vars_source
      else if (vars_source == "PUXURS" || vars_source == "PUXUYS") {
        for (int inode = 0; inode < num_nodes_Tecplot; inode++) {

          double entropy = var_source[inode][FOURTH];
          double pressure = var_source[inode][FIRST];
          double ux = var_source[inode][SECOND];
          double uy = var_source[inode][THIRD];

          var_source[inode][FIRST] = entropy;
          var_source[inode][SECOND] = ux;
          var_source[inode][THIRD] = uy;
          var_source[inode][FOURTH] = pressure;

        } // inode
      } // vars_source
      if (izone == num_zones_source - 1)
        mpi::wait_allothers("The solution to be interpolated is converted to be consistent with the current solution.");

    } // myinput->model_pde
    else
      mpi::graceful_exit("Unsupported physical model for the interpolation.");

  } // izone



  // inspect if every Tecplot element is right-handed
  for (int izone = 0; izone < num_zones_source; izone++) {

    // each zone may have different shape and number of control volumes
    int num_vertices_per_cv = TecplotFEZone[izone].num_vertices_per_cv;
    int num_cvs = TecplotFEZone[izone].num_cvs;
    ALLOCATE1D_INT_1ARG(nodeID_donor, num_vertices_per_cv);
    ALLOCATE2D_DOUBLE(xyz_donor, num_vertices_per_cv, num_dim_source);

    for (int icv = 0; icv < num_cvs; icv++) {

      // retrieve the node IDs of vertices constituting the current control volume, icv
      for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
        nodeID_donor[ivertex] = (TecplotFEZone[izone].connectivity[icv])[ivertex];

      // element handedness can be left-handed in some occasions
      // if so, a different permutation of vertex ordering is tried
      int handedness = NONE;
      do {

        // retrieve the donor vertices' xyz values
        for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
          for (int idim = XDIR; idim < num_dim_source; idim++)
            (xyz_donor[ivertex])[idim] = (TecplotFEZone[izone].xyz[nodeID_donor[ivertex]])[idim];

        for (int ivertex = 0; ivertex < num_vertices_per_cv - 1; ivertex++) {

          // get the indices of three neighboring vertices
          int ivertex0 = ivertex;
          int ivertex1 = (ivertex0 + 1) % num_vertices_per_cv;
          int ivertex2 = (ivertex0 + 2) % num_vertices_per_cv;

          // form two vectors
          for (int idim = XDIR; idim < num_dim_source; idim++) {

            vec0[idim] = (xyz_donor[ivertex1])[idim] - (xyz_donor[ivertex0])[idim];
            vec1[idim] = (xyz_donor[ivertex2])[idim] - (xyz_donor[ivertex1])[idim];

          } // idim

          // see if everything is right-handed
          if (math_matrix::cross_product(vec0, vec1, num_dim_source) < 0.0) {

            handedness = LEFT;
            std::cout << ">>> An element has a left-handedness for zone " << std::setw(3) << izone
                      << ", xyz = " << (xyz_donor[ivertex0])[XDIR] << ", "
                                    << (xyz_donor[ivertex0])[YDIR] << "; vertex order is reversed." << std::endl;
            break; // break the ivertex loop and move on to a new permutation in nodeID_donor

          } // math_matrix::cross_product(vec0, vec1, num_dim_source)
          else
            handedness = RIGHT;
        } // ivertex
        if (handedness == RIGHT) break; // if this element is right-handed, no need to check a further permutation

      } while ( std::next_permutation(nodeID_donor,nodeID_donor+num_vertices_per_cv) );

      // we are doomed..
      if (handedness == LEFT) {
        std::cout << "Left-handed element is not removed in zone " << std::setw(3) << izone << ", node ID ";
        for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
          std::cout << (TecplotFEZone[izone].connectivity[icv])[ivertex] << " ";
        std::cout << std::endl;
        mpi::graceful_exit("The right handedness of Tecplot element cannot be achieved; better to regenerate it.");
      } // handedness

      // in case handedness is fixed, override connectivity
      for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++) {
        (TecplotFEZone[izone].connectivity[icv])[ivertex] = nodeID_donor[ivertex];
      } // ivertex

    } // icv
    DEALLOCATE_1DPTR(nodeID_donor);
    DEALLOCATE_2DPTR(xyz_donor, num_vertices_per_cv);

  } // izone
  mpi::wait_allothers();
  mpi::wait_allothers("The element handedness is checked.");



  // loop over the entire interpolees (e.g. PLOT3D points) that the current core has and do interpolation
  if (myinput->num_dim != 2)
    mpi::graceful_exit("The following algorithm will not work for other than 2-D configurations.");
  double *matrixA = new double[num_dim_source * num_dim_source];
  double *vectorx = new double[num_dim_source];
  double *vectorb = new double[num_dim_source];
  for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
    for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
      for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
        int l0 = mygrid->idx1D(i, j, k);

        for (int idim = XDIR; idim < DIM_MAX; idim++)
          xyz_receiver[idim] = mygrid->cell[l0].xyz[idim];

        int found = FALSE; // by default, no cell is found for this point l0

        // start to search a matching control volume in every source zone
        for (int izone = 0; izone < num_zones_source; izone++) {

          // each zone may have different shape and number of control volumes
          int num_vertices_per_cv = TecplotFEZone[izone].num_vertices_per_cv;
          int num_cvs = TecplotFEZone[izone].num_cvs;
          int num_nodes = TecplotFEZone[izone].num_nodes;
          ALLOCATE1D_INT_1ARG(nodeID_donor, num_vertices_per_cv);
          ALLOCATE2D_DOUBLE(xyz_donor, num_vertices_per_cv, num_dim_source);
          ALLOCATE2D_DOUBLE(var_donor, num_vertices_per_cv, num_vars_source);

          // loop over the control volumes within this zone either until a 
          // corresponding control volume is found, or until control volumes 
          // are exhausted for search
          int icv = 0;
          while ( !found && icv < num_cvs ) {

            // retrieve the node IDs of vertices constituting the current control volume, icv
            for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
              nodeID_donor[ivertex] = (TecplotFEZone[izone].connectivity[icv])[ivertex];

            // retrieve the donor vertices' xyz values
            for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
              for (int idim = XDIR; idim < num_dim_source; idim++)
                (xyz_donor[ivertex])[idim] = (TecplotFEZone[izone].xyz[nodeID_donor[ivertex]])[idim];

            // check if this cell is where I am
            // case 1: if the current interpolation point lies exactly on top of 
            //         one of the verticies
            for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++) {

              for (int idim = XDIR; idim < num_dim_source; idim++)
                vec0[idim] = xyz_receiver[idim] - (xyz_donor[ivertex])[idim];
              if (sqrt(math_matrix::inner_product(vec0, vec0, num_dim_source)) <= eps_proximity)
                found = TRUE;

            } // ivertex

            // case 2:
            // a strategy here is to form a pair of vectors (both having the origins 
            // at a common vertex, but one having its destination at the neighboring vertex, 
            // while the other having its destination at the interpolation point)
            // then, a cross product is computed; if the signs for all the pairs are 
            // positive, the interpolation point is within this control volume
            // (since right-handedness is already checked above); otherwise, it is not; 
            // be careful about cases where the cross product is zero; the two vectors 
            // are co-linear then
            if ( !found ) {

              int point_inside = TRUE;
              int ivertex = 0;
              while ( point_inside && ivertex < num_vertices_per_cv ) {

                // get the indices of two neighboring vertices
                int ivertex0 = ivertex;
                int ivertex1 = (ivertex0 + 1) % num_vertices_per_cv;

                // form a pair of vectors mentioned above, vec0 and vec1
                for (int idim = XDIR; idim < num_dim_source; idim++) {

                  vec0[idim] = (xyz_donor[ivertex1])[idim] - (xyz_donor[ivertex0])[idim];
                  vec1[idim] = xyz_receiver[idim]          - (xyz_donor[ivertex0])[idim];

                } // idim
                double mag_vec0 = math_matrix::inner_product(vec0, vec0, num_dim_source);
                double mag_vec1 = math_matrix::inner_product(vec1, vec1, num_dim_source);

                // take a cross product (!!!this implementation is limited only to 2-D!!!)
                // note that this cross_product is equivalent to a sine of the angles between the two vectors
                double cross_product = math_matrix::cross_product(vec0, vec1, num_dim_source) / sqrt(mag_vec0 * mag_vec1);
                if (cross_product < -eps_smallAngle) // case 2-1: the interpolation point lies outside of this control volume
                  point_inside = FALSE;

//                else if (abs(cross_product) <= eps_smallAngle) { // case 2-2: the two vectors are co-linear
//
//                  for (int idim = XDIR; idim < num_dim_source; idim++)
//                    vec0[idim] = xyz_receiver[idim] - (xyz_donor[ivertex1])[idim];
//                  double inner_product = math_matrix::inner_product(vec0, vec1, num_dim_source);
//
//                  if (inner_product > 0.0)
//                    point_inside = FALSE; // co-linear but outside of this control volume
//                                          // if inner_product < 0, the point is between the two vertices
//                  else
//                    point_inside = TRUE;
//
//                  ivertex = num_vertices_per_cv; // enforce the exit condition
//
//                } // abs(cross_product)

                else if (cross_product >= -eps_smallAngle) // case 2-3: the point is still inside; keep searching
                  point_inside = TRUE;

                ivertex++;

              } // point_inside

              if ( point_inside ) // a control volume enclosing our recipient is found
                found = TRUE;

              icv++;

            } // !found

            if ( found ) {

              // retrieve the donor nodes' variables to be interpolated
              for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
                for (int ivar = 0; ivar < num_vars_source; ivar++)
                  (var_donor[ivertex])[ivar] = (TecplotFEZone[izone].var[nodeID_donor[ivertex]])[ivar];

              // the following linear model fails near x = 0, y = 0, z = 0 since the coefficient matrix gets singular
              // a quickest fix is to translate the receiver and donor locations by the same amount
              for (int idim = XDIR; idim < num_dim_source; idim++)
                xyz_receiver[idim] -= xyz_min[idim] - 1.0; // set the minimum location arbitrarily to 1.0 in every direction
              for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
                for (int idim = XDIR; idim < num_dim_source; idim++)
                  (xyz_donor[ivertex])[idim] -= xyz_min[idim] - 1.0; // set the minimum location arbitrarily to 1.0 in every direction

              // construct the coefficient matrix in 2-D
              // this is done just once for this interpolation point since the elements are purely geometrical
              for (int irow = XDIR; irow < num_dim_source; irow++)
                for (int icol = XDIR; icol < num_dim_source; icol++)
                  matrixA[ij2by2[irow][icol]] = 0.0;
              for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++) {

                matrixA[ij2by2[0][0]] += (xyz_donor[ivertex])[XDIR] * (xyz_donor[ivertex])[XDIR];
                matrixA[ij2by2[0][1]] += (xyz_donor[ivertex])[XDIR] * (xyz_donor[ivertex])[YDIR];
                matrixA[ij2by2[1][1]] += (xyz_donor[ivertex])[YDIR] * (xyz_donor[ivertex])[YDIR];

              } // ivertex
              matrixA[ij2by2[1][0]] = matrixA[ij2by2[0][1]]; // enforce symmetry

              for (int ivar = 0; ivar < num_vars_source; ivar++) {

                // right-hand-side vector is dependent on variable
                for (int icol = XDIR; icol < num_dim_source; icol++)
                  vectorb[icol] = 0.0;
                for (int ivertex = 0; ivertex < num_vertices_per_cv; ivertex++)
                  for (int icol = XDIR; icol < num_dim_source; icol++)
                    vectorb[icol] += (xyz_donor[ivertex])[icol] * (var_donor[ivertex])[ivar];

                // solve the n-by-n linear system for the coefficients
                double determinant = math_matrix::solve_linearSystem(matrixA, vectorx, vectorb, num_dim_source);

                // take a linear combination
                (mystate->sol[ivar])[l0] = 0.0;
                for (int idim = XDIR; idim < num_dim_source; idim++)
                  (mystate->sol[ivar])[l0] += vectorx[idim] * xyz_receiver[idim];

              } // ivar
            } // found
          } // !found

          DEALLOCATE_1DPTR(nodeID_donor);
          DEALLOCATE_2DPTR(xyz_donor, num_vertices_per_cv);
          DEALLOCATE_2DPTR(var_donor, num_vertices_per_cv);

        } // izone

        // in case that no enclosing control volume is found
        if ( !found ) {

          for (int ivar = FIRST; ivar < num_vars_source; ivar++)
            (mystate->sol[ivar])[l0] = DUMMY_DOUBLE;
          std::cout << "x = " << xyz_receiver[XDIR] << ", y = " << xyz_receiver[YDIR] << " doesn't have any enclosing cell at core " << mpi::irank << ", l0: " << l0 << "; go with the nearest neighbor." << std::endl;

          // do a nearest-neighbor interpolation
          // find the minimum distance to the donor points by a brute force
          int izone_distance_min = NONE;
          int inode_distance_min = NONE;
          double distance_min = DUMMY_LARGE;
          for (int izone = 0; izone < num_zones_source; izone++) {

            // each zone may have different shape and number of control volumes
            int num_nodes = TecplotFEZone[izone].num_nodes;
            double **xyz_source = TecplotFEZone[izone].xyz;
            double **var_source = TecplotFEZone[izone].var;

            for (int inode = 0; inode < num_nodes; inode++) {

              double distance = 0.0;
              for (int idim = XDIR; idim < num_dim_source; idim++)
                distance += pow(xyz_receiver[idim] - (xyz_source[inode])[idim], 2); // sqrt is not applied since only the closest-point's index is of interest

              if (distance < distance_min) {

                izone_distance_min = izone;
                inode_distance_min = inode;
                distance_min = distance;

              } // distance
            } // inode
          } // izone
          assert( izone_distance_min != NONE );
          assert( inode_distance_min != NONE );
          assert( distance_min != DUMMY_LARGE );
//          distance_min = sqrt( distance_min ); // the minimum distance itself is not used

          // get the solution of the nearest neighbor
          for (int ivar = 0; ivar < num_vars_source; ivar++)
            (mystate->sol[ivar])[l0] = (TecplotFEZone[izone_distance_min].var[inode_distance_min])[ivar];

        } // !found

      } // i
    } // j
  } // k
  DEALLOCATE_1DPTR(matrixA);
  DEALLOCATE_1DPTR(vectorx);
  DEALLOCATE_1DPTR(vectorb);
  std::cout << ">>> Rank " << std::setw(6) << mpi::irank << " out of " << std::setw(6) << mpi::nprocs << " is done with interpolation." << std::endl;
  mpi::wait_allothers();
  mpi::wait_allothers("The linear interpolation is completed.");



  // adjust the solution-vector ordering for 2-D cases
  if (myinput->num_dim == 2) {

    for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
      for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
        for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
          int l0 = mygrid->idx1D(i, j, k);

          double tmp = (mystate->sol[IVAR_P])[l0];
          (mystate->sol[IVAR_P])[l0] = (mystate->sol[IVAR_UZ])[l0];
          (mystate->sol[IVAR_UZ])[l0] = tmp;

        } // i
      } // j
    } // k
    mpi::wait_allothers("Since it's a 2-D case, the order of the solution vector is adjusted.");

  } // myinput->num_dim



  // apply a filter
  for (int ifilter = 0; ifilter < myinput->num_filters_interpolation; ifilter++) {

    spatial::apply_filter(mystate->num_vars_sol, mystate->num_samples, mystate->sol, mygrid);
    if (mpi::irank == 0)
      if (ifilter == 0)
        std::cout << ">>> The filter is applied " << std::setw(6) << ifilter + 1 << " time." << std::endl;
      else
        std::cout << ">>> The filter is applied " << std::setw(6) << ifilter + 1 << " times." << std::endl;

  } // ifilter



  // write the result
  mystate->time_step = 0;
  mystate->time_sol = 0.0;
  if (myinput->interpolate_into_which_PLOT3D == "PLOT3D_SOLUTION") {

    io::write_solution(myinput, mygrid, block, mystate);
    mpi::wait_allothers("The interpolated solution is written; check the regular solution file written at the counter 0.");

  } // myinput->interpolate_into_which_PLOT3D
  else if (myinput->interpolate_into_which_PLOT3D == "PLOT3D_FUNCTION") {

    for (int k = mygrid->iso[ZETA]; k <= mygrid->ieo[ZETA]; k++) {
      for (int j = mygrid->iso[ETA]; j <= mygrid->ieo[ETA]; j++) {
        for (int i = mygrid->iso[XI]; i <= mygrid->ieo[XI]; i++) {
          int l0 = mygrid->idx1D(i, j, k);

          for (int ivar = FIRST; ivar < mystate->num_vars_sol; ivar++)
            (mystate->sol_mean[ivar])[l0] = (mystate->sol[ivar])[l0];

        } // i
      } // j
    } // k
    mpi::wait_allothers();

    io::write_solution_mean(myinput, mygrid, block, mystate);
    mpi::wait_allothers("The interpolated mean state is written; check the solution file starting with 'mean_'.");

  } // myinput->interpolate_into_which_PLOT3D

  mpi::end_parallel_global();

  return 0;

} // main
