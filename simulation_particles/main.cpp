#include "main_solver.h"

int main(){
  int Nx = 1; //number of layers in x
  int max_N_ = 10;

  for(int n = 1; n <= max_N_; n += 1){
    simulate_n_particles simulator(n, Nx);
    simulator.main();
  }

  return 0;
}
