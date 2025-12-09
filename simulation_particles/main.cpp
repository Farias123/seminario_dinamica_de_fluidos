#include "main_solver.h"

int main(int argc, char *argv[]){
  int N = stoi(argv[1]);
  bool save_positions = true;

  simulate_n_particles simulator(N, save_positions);
  simulator.main();

  return 0;
}
