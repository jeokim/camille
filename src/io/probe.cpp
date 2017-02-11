#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include "probe.h"

namespace probe {

void initialize(UserInput *myinput) {

  MESSAGE_STDOUT("Hi, my name is Youngchul Kim.");
  mpi::graceful_exit("See you later!");

  return;

} // initialize

} // probe
