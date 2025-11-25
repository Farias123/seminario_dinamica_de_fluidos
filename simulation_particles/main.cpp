#include "main_solver.h"

int main(){
  int max_N = 10;

  for(int n = 5; n <= max_N; n += 1){
    simulate_n_particles simulator(n);
    simulator.main();
  }

  return 0;
}
