#include "main_solver.h"

int main(){
  int max_N = 15000;

  for(int n = 15.000; n <= max_N; n += 1){
    simulate_n_particles simulator(n);
    simulator.main();
  }

  return 0;
}
