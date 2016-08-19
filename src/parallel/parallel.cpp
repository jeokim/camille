#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "../core/param.h"
#include "../core/macros_inlines.h"
#include "../math/constants.h"
#include "../math/algebra.h"
#include "parallel.h"

namespace mpi {

int irank;
int nprocs;

int ierr;

MPI_Status *status = MPI_STATUS_IGNORE;
MPI_Request request[2];

// for region-level message passing
MPI_Comm comm_region;
MPI_Group group_region;
int irank_region;
int nprocs_region;

// for block-level message passing
MPI_Comm comm_block;
MPI_Group group_block;
int irank_block;
int nprocs_block;

// for message passing amount head ranks
MPI_Comm comm_head;
MPI_Group group_head;
int irank_head;
int nprocs_head;

void wait_allothers() {

  MPI_Barrier(MPI_COMM_WORLD);

} // wait_allothers



void wait_allothers(std::string wait_message) {

  MPI_Barrier(MPI_COMM_WORLD);
//  sleep(0.5);
  MESSAGE_STDOUT(wait_message);

} // wait_allothers



void wait_cores_in_(MPI_Comm comm_in) {

  MPI_Barrier(comm_in);

} // wait_cores_in_



void init_parallel_global() {

  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  wait_allothers("Global MPI initiated.");

} // init_parallel_global



void end_parallel_global() {

  wait_allothers("Global MPI terminated.");

//  MPI_Abort(MPI_COMM_WORLD, ierr);
  MPI_Finalize();
  std::exit(0);

} // end_parallel_global



void graceful_exit(std::string exit_message) {

  wait_allothers(exit_message);

  end_parallel_global();

} // graceful_exit

} // parallel
