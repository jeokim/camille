#include "camille.h"

int main(int argc, char * argv[]) {

  simulation::initialize(argc, argv);
  simulation::run();
  simulation::finalize();

} // main
