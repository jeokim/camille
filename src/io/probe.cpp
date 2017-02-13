#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "probe.h"

namespace probe {

int num_myprobes;
ProbePoint *probe_point;

void initialize(UserInput *myinput, Geometry::StructuredGrid *mygrid) {

  if (myinput->do_probe == FALSE)
    return;
  assert(myinput->num_probes > 0);

  int *corresponding_ijk;
  std::vector<int> ijk_local;
  std::stringstream str_dummy;
  std::string str_output;

  num_myprobes = 0; // no probe for now

  // core2probe[k] = 1 if this core owns the probe k (k = 0, 1, 2, ...)
  int *core2probe;
  ALLOCATE1D_INT_1ARG(core2probe, myinput->num_probes);

  // go over all probes
  ALLOCATE1D_INT_1ARG(corresponding_ijk, DIM_MAX);
  for (int iprobe = 0; iprobe < myinput->num_probes; iprobe++) {

    // x, y, & z locations of a point probe
    double xyz[DIM_MAX];
    for (int idir = XDIR; idir < DIM_MAX; idir++)
      xyz[idir] = myinput->tmp_probe_xyz[iprobe][idir];

    // search
    if (mygrid->check_if_this_is_my_point(myinput->num_dim, xyz, corresponding_ijk) == TRUE) {

      num_myprobes++;
      core2probe[iprobe] = TRUE; // the current core owns this probe (for now)
      for (int idir = XI; idir < DIM_MAX; idir++)
        ijk_local.push_back(corresponding_ijk[idir]); // keep the corresponding cell's grid-level ijk indices (for now)

    } // mygrid->check_if_this_is_my_point(myinput->num_dim, xyz, corresponding_ijk)
  } // iprobe
  DEALLOCATE_1DPTR(corresponding_ijk);

  // due to grid overlapping (either ghost cell or overset), a single probe could be claimed by more than one grid
for (int iprobe = 0; iprobe < myinput->num_probes; iprobe++)
std::cout << "Rank: " << mpi::irank << ", probe: " << iprobe << ", num probes: " << num_myprobes << ", core2probe: " << core2probe[iprobe] << std::endl;
mpi::graceful_exit("bye!");
  int *core2probe_sum;
  ALLOCATE1D_INT_1ARG(core2probe_sum, myinput->num_probes);
  MPI_Allreduce(core2probe, core2probe_sum, myinput->num_probes, MPI_INT, MPI_SUM, mpi::comm_region);

  for (int iprobe = 0; iprobe < myinput->num_probes; iprobe++) {

    // case 1: no core claims this probe
    if (core2probe_sum[iprobe] == 0) {

      for (int idir = XDIR; idir < DIM_MAX; idir++) {
        str_dummy.str("");
        str_dummy << myinput->tmp_probe_xyz[iprobe][idir];
        str_output += str_dummy.str() + " ";
      } // idir
      mpi::graceful_exit("Probe (x y z) = ( " + str_output + ") does not belong to any core.");

    } // core2probe_sum[iprobe]
    // case 2: there is a single core which owns this probe
    else if (core2probe_sum[iprobe] == 1) {

      // nothing to do

    } // core2probe_sum[iprobe]
    // case 3: there are more than a single core which owns this probe
    else if (core2probe_sum[iprobe] > 1) {

      int claimed = FALSE; // whether this core has claimed this probe or not
      int claimed_sum = 0;

      // pick the one owned by the lowest rank
      for (int irank = 0; irank < mpi::nprocs; irank++) {
        if (mpi::irank == irank) {

          if (core2probe[iprobe] == TRUE) {
            claimed = TRUE;
          } // core2probe[iprobe]

        } // mpi::irank
        mpi::wait_allothers();
        MPI_Allreduce(&claimed, &claimed_sum, 1, MPI_INT, MPI_SUM, mpi::comm_region);
        if (claimed_sum == TRUE) { // if some core (including mine) has already taken this probe
          if (claimed == FALSE) { // but my core has not; thus, give up this probe
            num_myprobes--;
            core2probe[iprobe] = FALSE;
            for (int idir = XI; idir < DIM_MAX; idir++)
              ijk_local[iprobe*DIM_MAX+idir] = NONE;
          } // claimed
        } // claimed_sum
      } // irank
    } // core2probe_sum[iprobe]
  } // iprobe

  // ensure that a probe is taken a single core and only by a single core
  MPI_Allreduce(core2probe, core2probe_sum, myinput->num_probes, MPI_INT, MPI_SUM, mpi::comm_region);
  for (int iprobe = 0; iprobe < myinput->num_probes; iprobe++) {
    if (core2probe_sum[iprobe] != 1) {
      for (int idir = XDIR; idir < DIM_MAX; idir++) {
        str_dummy << myinput->tmp_probe_xyz[iprobe][idir];
        str_output += str_dummy.str()+" ";
      } // idir
      mpi::graceful_exit("Probe x,y,z = " + str_output + "is not properly store; check if everything is okay.");
    } // core2probe_sum[iprobe]
  } // iprobe

  // store
  if (num_myprobes > 0) {
    probe_point = new ProbePoint[num_myprobes];

    int counter = 0;
    for (int iprobe = 0; iprobe < myinput->num_probes; iprobe++) {
      if (core2probe[iprobe] == TRUE) {

        probe_point[counter].name = myinput->tmp_probe_name[iprobe];
        probe_point[counter].interval = myinput->tmp_probe_interval[iprobe];
        for (int idir = XDIR; idir < DIM_MAX; idir++)
          probe_point[counter].xyz[idir] = myinput->tmp_probe_xyz[iprobe][idir];
        for (int idir = XI; idir < DIM_MAX; idir++)
          probe_point[counter].ijk[idir] = ijk_local[iprobe*DIM_MAX+idir];
        ALLOCATE1D_DOUBLE_1ARG(probe_point[counter].fac_interp, static_cast<int>(pow(2,myinput->num_dim)));

        // precompute interpolation factors




        counter++;

      } // core2probe[iprobe]
    } // iprobe
    assert(counter == num_myprobes);
  } // num_myprobes

  // clean up temporary storage
  ijk_local.clear();
  DEALLOCATE_1DPTR(core2probe);
  DEALLOCATE_1DPTR(core2probe_sum);
  DEALLOCATE_1DPTR(myinput->tmp_probe_name);
  DEALLOCATE_1DPTR(myinput->tmp_probe_interval);
  DEALLOCATE_2DPTR(myinput->tmp_probe_xyz, myinput->num_probes);

mpi::graceful_exit("Probe-module test completed.\n");

  return;

} // initialize

} // probe
