#include "main_solver.h"

int main(){
  int Nx = 10; //number of layers in x
  int max_N_ = 3;

  for(int n = 3; n <= max_N_; n += 1){
    simulate_n_particles simulator(n, Nx);
    simulator.main();
  }

  return 0;
}
