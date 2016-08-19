#ifndef PARALLEL_PARALLEL_H
#define PARALLEL_PARALLEL_H

//

#include <string>
#include <mpi.h>

namespace mpi {

extern int irank;
extern int nprocs;

extern int ierr;

extern MPI_Status *status;
extern MPI_Request request[2];

// for region-level message passing
extern MPI_Comm comm_region; // equivalent to MPI_COMM_WORLD
extern MPI_Group group_region;
extern int irank_region; // equivalent to irank
extern int nprocs_region; // equivalent to nprocs

// for block-level message passing
extern MPI_Comm comm_block;
extern MPI_Group group_block;
extern int irank_block;
extern int nprocs_block;

// for message passing amount head ranks
extern MPI_Comm comm_head;
extern MPI_Group group_head;
extern int irank_head;
extern int nprocs_head;

void wait_allothers(void);
void wait_allothers(std::string);
void wait_cores_in_(MPI_Comm);

void init_parallel_global(void);
void end_parallel_global(void);
void graceful_exit(std::string);

} // mpi

//

#endif
