#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "probe.h"

namespace probe {

void initialize(UserInput *myinput) {

  if (myinput->do_probe == FALSE)
    return;
  assert(myinput->num_probes > 0);











//     for (int iprobe = 0; iprobe < num_probes; iprobe++) {
//       inputDeck::get_userInput("PROBE","NAME",tmp_probe_name[iprobe],iprobe);
//       inputDeck::get_userInput("PROBE","INTERVAL",tmp_probe_interval[iprobe],iprobe);
//       inputDeck::get_userInput("PROBE","XYZ",DIM_MAX,xyz,iprobe);
//       for (int idir = XDIR; idir < DIM_MAX; idir++)
//         tmp_probe_xyz[iprobe][idir] = xyz[idir];
//     } // iprobe

  // clean up temporary storage  
  DEALLOCATE_1DPTR(myinput->tmp_probe_name);
  DEALLOCATE_1DPTR(myinput->tmp_probe_interval);
  DEALLOCATE_2DPTR(myinput->tmp_probe_xyz, myinput->num_probes);

  return;

} // initialize

} // probe
