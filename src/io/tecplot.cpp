#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "tecplot.h"

namespace tecplot {

const int line_width = 200;

void read_ASCII_header_FE(FILE *pFile, int num_dim, int num_vars) {

  char str[line_width];

  // TITLE
  fgets(str, line_width, pFile);
//  if (mpi::irank == 0)
//    std::cout << str;

  // VARIABLES
  for (int idim = XDIR; idim < num_dim; idim++) {

    fgets(str, line_width, pFile);
//    if (mpi::irank == 0)
//      std::cout << str;

  } // idim
  for (int ivar = FIRST; ivar < num_vars; ivar++) {

    fgets(str, line_width, pFile);
//    if (mpi::irank == 0)
//      std::cout << str;

  } // ivar

  return;

} // read_ASCII_header_FE



int read_ASCII_zoneHeader_FE(FILE *pFile, int num_dim, int num_vars, int &num_nodes, int &num_cvs) {

  char str[line_width];
  std::string str_tmp;
  int num_vertices_per_cv = 0;

  // ZONE
  fgets(str, line_width, pFile);
//  if (mpi::irank == 0)
//    std::cout << str;

  fgets(str, line_width, pFile);
//  if (mpi::irank == 0)
//    std::cout << str;

  fscanf(pFile, " Nodes=%d, Elements=%d, ZONETYPE=%s", &num_nodes, &num_cvs, str);

  str_tmp = std::string(str);
  if (str_tmp == "FETriangle")
    num_vertices_per_cv = 3;
  else if (str_tmp == "FEQuadrilateral")
    num_vertices_per_cv = 4;
  else
    mpi::graceful_exit("Only FETriangle and/or FEQuadrilateral format is supported for the Tecplot ASCII file-reading; you have \n" + str_tmp);

  if (mpi::irank == 0)
    std::cout << ">>> Number of nodes: "  << std::setw(8) << num_nodes
              << ", number of elements: " << std::setw(8) << num_cvs
              << ", type of elements: "   << str_tmp << std::endl;

  // DATAPACKING
  fscanf(pFile, " DATAPACKING=%s", str);
  str_tmp = std::string(str);
  if (str_tmp != "POINT")
    mpi::graceful_exit("Only POINT format is supported for the Tecplot ASCII file-reading.");
  if (mpi::irank == 0)
    std::cout << ">>> Data packing is " << str_tmp << std::endl;

  // somehow I need this additional fgets, presumably to read a carriage return not read by fscanf?
  fgets(str, line_width, pFile);
//  if (mpi::irank == 0)
//    std::cout << str;

  fgets(str, line_width, pFile);
//  if (mpi::irank == 0)
//    std::cout << str;

  return num_vertices_per_cv;

} // read_ASCII_zoneHeader_FE



void read_ASCII_POINT(FILE *pFile, int num_dim, int num_vars, int num_nodes, int num_cvs, double ** &xyz, double ** &var, int ** &connectivity, int num_vertices_per_cv) {

  char str[line_width];

  // first allocate
  ALLOCATE2D_DOUBLE(xyz, num_nodes, num_dim);
  ALLOCATE2D_DOUBLE(var, num_nodes, num_vars);
  ALLOCATE2D_INT(connectivity, num_cvs, num_vertices_per_cv);

  // read the nodal variables
  // for some reasons, fscanf appears not to work well in the current parallel environment
  // since the Tecplot file is read just once, let all the cores read the file sequentially
  if (mpi::irank == 0)
    std::cout << ">>> Nodal data are being read: ";
  for (int irank = FIRST; irank < mpi::nprocs; irank++) {
    if (mpi::irank == irank) {

      // field variables
      for (int inode = FIRST; inode < num_nodes; inode++) {

        for (int idim = XDIR; idim < num_dim; idim++)
          fscanf(pFile, "%lf", &(xyz[inode][idim]));

        for (int ivar = FIRST; ivar < num_vars; ivar++)
          fscanf(pFile, "%lf", &(var[inode][ivar]));

      } // inode

      // connectivity for control volume
      for (int icv = FIRST; icv < num_cvs; icv++)
        for (int ivertex = FIRST; ivertex < num_vertices_per_cv; ivertex++)
          fscanf(pFile, "%d", &(connectivity[icv][ivertex]));

      std::cout << ".";
    } // mpi::irank
    mpi::wait_allothers();

  } // irank
  if (mpi::irank == 0)
    std::cout << std::endl;
//  if (mpi::irank == 0) {
//    for (int inode = FIRST; inode < num_nodes; inode++) {
//
//      for (int ivar = FIRST; ivar < num_dim; ivar++)
//        std::cout << " " << std::setw(20) << std::scientific << std::setprecision(18) << xyz[inode][ivar];
//
//      for (int ivar = FIRST; ivar < num_vars; ivar++) {
//        std::cout << " " << std::setw(20) << std::scientific << std::setprecision(18) << var[inode][ivar];
//        if (ivar == num_vars - 1)
//          std::cout << std::endl;
//      } // ivar
//    } // inode
//    std::cout << "x: " << xyz[139936-1][XDIR] << "; y: " << xyz[139936-1][YDIR] << std::endl;
//    std::cout << "vars: " << var[139936-1][0] << "; " << var[139936-1][1] << "; " << var[139936-1][2] << "; " << var[139936-1][3] << std::endl;
//    for (int icv = FIRST; icv < num_cvs; icv++)
//      for (int inode = FIRST; inode < num_vertices_per_cv; inode++) {
//        std::cout << " " << connectivity[icv][inode];
//        if (inode == num_vertices_per_cv - 1)
//          std::cout << std::endl;
//      } // inode
//    std::cout << "connectivity: " << connectivity[278153-1][0] << " " << connectivity[278153-1][1] << " " << connectivity[278153-1][2] << std::endl;
//  } // mpi::irank

  // since this is a C++ code, decrease what's in the connectivity since indices start from zero
  for (int icv = FIRST; icv < num_cvs; icv++)
    for (int ivertex = FIRST; ivertex < num_vertices_per_cv; ivertex++)
      connectivity[icv][ivertex]--;

  // report the minimum and maximum values of each variable
  for (int idim = XDIR; idim < num_dim; idim++) {
    double min = DUMMY_LARGE, max = -DUMMY_LARGE;

    for (int inode = FIRST; inode < num_nodes; inode++) {

      min = std::min(min, xyz[inode][idim]);
      max = std::max(max, xyz[inode][idim]);

    } // inode
    if (mpi::irank == 0)
      std::cout << ">>> min/max of variable " << std::setw(3) << idim << " of xyz are " << min << " / " << max << std::endl;

  } // idim
  //
  for (int ivar = FIRST; ivar < num_vars; ivar++) {
    double min = DUMMY_LARGE, max = -DUMMY_LARGE;

    for (int inode = FIRST; inode < num_nodes; inode++) {

      min = std::min(min, var[inode][ivar]);
      max = std::max(max, var[inode][ivar]);

    } // inode
    if (mpi::irank == 0)
      std::cout << ">>> min/max of the solution variable " << std::setw(3) << ivar << " are " << min << " / " << max << std::endl;

  } // ivar

  // somehow I need this additional fgets, presumably to read a carriage return not read by fscanf?
  fgets(str, line_width, pFile);
//  if (mpi::irank == 0)
//    std::cout << str;

  return;

} // read_ASCII_POINT

} // tecplot
