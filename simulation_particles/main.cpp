#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

const double dt = 1.0/30;
const double initial_speed_u = 1.0; // m/s

//atom variables
const double n_avogadro = 6.0221408e23;
const double mass_he = 0.004002602/n_avogadro; //kg
const double radius_he = 1.4e-10; //m
const double density_he = 0.1785; //kg/m3

//space between atoms surfaces and their square limit
const double r_ = pow(mass_he/(8*density_he), 1.0/3.0) - radius_he;


void update_positions(){
// create class to solve verlet
  cout << "Verlet + colisões";
  //   save_position_file();
}

void save_position_file(){
  cout << "Baseado no número de partículas simuladas, se for menos que 10 mil, salva os dados de todas. Acima disso"
  << "sorteia 10 mil números aleatórios entre 0 e Nx*N_*N_. Se o índice da partícula (i + Nx*j + Nx*N_*k) é um dos"
  << "sorteados, será salvo no arquivo do passo de tempo. Esquema: criar uma pasta com nome tendo os N_ e Nx, para cada"
  << "passo criar arquivo com posições de todas as partículas";
}


double simulate_n_particles(int N_, int Nx){
  //N = Nx*N_^2

  const double sphere_radius = N_*(radius_he + r_)/2.0;

  auto idx = [N_, Nx](int i, int j, int k, int dim){
    return  (N_*(N_*i + j) + k)*3 + dim;
  };

  cout << "N_ = "<< N_<<"\nR = " << sphere_radius << "\n";

  const double x0 = radius_he + r_ - 4.0*sphere_radius, y0 = radius_he + r_ - 2.0*sphere_radius, z0 = radius_he + r_ - 2.0*sphere_radius;

  double* positions = new double[Nx*N_*N_*3];
  double* velocities = new double[Nx*N_*N_*3];

  // initial conditions
  for(int i = 0; i < Nx; i += 1){
    for(int j = 0; j < N_; j += 1){
      for(int k = 0; k < N_; k += 1){
        positions[idx(i, j, k, 0)] = x0 - 2*(radius_he + r_)*i;
        positions[idx(i, j, k, 1)] = y0 + 2*(radius_he + r_)*j;
        positions[idx(i, j, k, 2)] = z0 + 2*(radius_he + r_)*k;

        velocities[idx(i, j, k, 0)] = initial_speed_u;
        velocities[idx(i, j, k, 1)] = 0.0;
        velocities[idx(i, j, k, 2)] = 0.0;
//        save_position_file();
      }
    }
  }

//  update_positions(positions, velocities);

  delete[] positions;
  delete[] velocities;

  return 0.0;
}

int main(){
  int Nx = 1; //number of layers in x
  int max_N_ = 10;

  for(int n = 1; n <= max_N_; n += 1){
    simulate_n_particles(n, Nx);
  }

  return 0;
}
