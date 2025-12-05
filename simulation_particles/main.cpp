#include "main_solver.h"

int main(){
  int max_N = 10000;

  for(int n = 10000; n <= max_N; n += 1){
    simulate_n_particles simulator(n);
    simulator.main();
  }

  return 0;
}
