#include "main_solver.h"

int main(){
  int max_N = 10000;
  bool save_positions = false;

  for(int n = 10000; n <= max_N; n += 1){
    simulate_n_particles simulator(n, save_positions);
    simulator.main();
  }

  return 0;
}
