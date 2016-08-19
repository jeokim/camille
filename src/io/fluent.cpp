#include <iostream>
#include <iomanip>
#include <fstream>

#include "fluent.h"

namespace fluent {

void read_profile(std::string file_profile_FLUENT, int num_dim, int num_vars, int &num_cvs, double ** &xyz, double ** &var) {

  std::ifstream ifs;
  std::string line_cur;

  ifs.open(cstr_to_constchar(file_profile_FLUENT), std::ifstream::in);
  if (!ifs.is_open())
    mpi::graceful_exit("The FLUENT profile to interpolate does not exist.");

  // step 1: count the number of control volumes in the profile
  std::getline(ifs, line_cur); // read the very 1st line of the file
  std::getline(ifs, line_cur); // read a variable label
  std::getline(ifs, line_cur); // read the 1st line of actual data

  num_cvs = 0;
  while (line_cur.at(FIRST) != ')') { // read until the right delimiter shows up

    num_cvs++;

    // read the next line whether the current line is skipped or not
    std::getline(ifs, line_cur);

  } // line_cur.at(FIRST)
//  std::cout << "Number of control volumes in the FLUENT calculation is " << std::setw(8) << num_cvs << std::endl;

  // step 2: allocate the arrays
  ALLOCATE2D_DOUBLE(xyz, num_cvs, num_dim);
  ALLOCATE2D_DOUBLE(var, num_cvs, num_vars);

  // step 3: read and store the pointwise data
  ifs.clear();
  ifs.seekg(0, ifs.beg);

  std::getline(ifs, line_cur); // read the very 1st line of the file

  // xyz first
  for (int idim = XDIR; idim < num_dim; idim++) {

    std::getline(ifs, line_cur); // read a variable label
    std::getline(ifs, line_cur); // read the 1st line of actual data

    int counter = 0;
    while (line_cur.at(FIRST) != ')') { // read until the corresponding right delimiter shows up

      xyz[counter][idim] = atof(cstr_to_constchar(line_cur));
//      if (mpi::irank == 0) {
//
//        std::cout << "Variable: " << std::setw(3) << idim
//                  << ", counter: " << std::setw(8) << counter
//                  << ", value: " << std::setprecision(10) << xyz[counter][idim] << std::endl;
//
//      } // mpi::irank

      counter++;

      // read the next line whether the current line is skipped or not
      std::getline(ifs, line_cur);

    } // line_cur.at(FIRST)
    assert( counter == num_cvs );

  } // idim

  // solution variables
  for (int ivar = 0; ivar < num_vars; ivar++) {

    std::getline(ifs, line_cur); // read a variable label
    std::getline(ifs, line_cur); // read the 1st line of actual data

    int counter = 0;
    while (line_cur.at(FIRST) != ')') { // read until the corresponding right delimiter shows up

      var[counter][ivar] = atof(cstr_to_constchar(line_cur));

      counter++;

      // read the next line whether the current line is skipped or not
      std::getline(ifs, line_cur);

    } // line_cur.at(FIRST)
    assert( counter == num_cvs );

  } // ivar
  ifs.close();

  // report the number of control volumes
  if (mpi::irank == 0)
    std::cout << ">>> The FLUENT profile contains " << num_cvs << " control volumes." << std::endl;

  // report the minimum and maximum values of each variable
  for (int idim = XDIR; idim < num_dim; idim++) {
    double min = DUMMY_LARGE, max = -DUMMY_LARGE;

    for (int icv = 0; icv < num_cvs; icv++) {

      min = std::min(min, xyz[icv][idim]);
      max = std::max(max, xyz[icv][idim]);

    } // icv
    if (mpi::irank == 0)
      std::cout << ">>> min/max of variable " << std::setw(3) << idim << " of xyz are " << min << " / " << max << std::endl;

  } // idim
  //
  for (int ivar = FIRST; ivar < num_vars; ivar++) {
    double min = DUMMY_LARGE, max = -DUMMY_LARGE;

    for (int icv = 0; icv < num_cvs; icv++) {

      min = std::min(min, var[icv][ivar]);
      max = std::max(max, var[icv][ivar]);

    } // icv
    if (mpi::irank == 0)
      std::cout << ">>> min/max of the solution variable " << std::setw(3) << ivar << " are " << min << " / " << max << std::endl;

  } // ivar

  return;

} // read_profile

} // fluent
