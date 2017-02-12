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

  // count how many probes this core contains
  num_myprobes = 0; // no probe for now

  // core2probe[k] = 1 if this core owns the probe k
  int *core2probe;
  ALLOCATE1D_INT_1ARG(core2probe,myinput->num_probes);

  // go over all probes
  for (int iprobe = 0; iprobe < myinput->num_probes; iprobe++) {

    // x, y, & z locations of a point probe
    double xyz[DIM_MAX];
    for (int idir = XDIR; idir < DIM_MAX; idir++)
      xyz[idir] = myinput->tmp_probe_xyz[iprobe][idir];

    // search
    for (int k = mygrid->is[ZETA]; k <= mygrid->ie[ZETA]; k++)
      for (int j = mygrid->is[ETA]; j <= mygrid->ie[ETA]; j++)
        for (int i = mygrid->is[XI]; i <= mygrid->ie[XI]; i++) {
  
          int l0 = mygrid->idx1D(i, j, k);



        } // i

  } // iprobe

  // store
  if (num_myprobes > 0) {
    probe_point = new ProbePoint[num_myprobes];

  } // num_myprobes

//     for (int iprobe = 0; iprobe < num_probes; iprobe++) {
//       inputDeck::get_userInput("PROBE","NAME",tmp_probe_name[iprobe],iprobe);
//       inputDeck::get_userInput("PROBE","INTERVAL",tmp_probe_interval[iprobe],iprobe);
//       inputDeck::get_userInput("PROBE","XYZ",DIM_MAX,xyz,iprobe);
//       for (int idir = XDIR; idir < DIM_MAX; idir++)
//         tmp_probe_xyz[iprobe][idir] = xyz[idir];
//     } // iprobe

  // clean up temporary storage
  DEALLOCATE_1DPTR(core2probe);
  DEALLOCATE_1DPTR(myinput->tmp_probe_name);
  DEALLOCATE_1DPTR(myinput->tmp_probe_interval);
  DEALLOCATE_2DPTR(myinput->tmp_probe_xyz,myinput->num_probes);

  return;

} // initialize

} // probe
