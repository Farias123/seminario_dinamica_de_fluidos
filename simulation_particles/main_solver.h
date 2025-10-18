#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

const double dt = 1e-10;
const double t_max = 1.0; //s

const double initial_speed_u = 1.0; // m/s

//atom variables
const double n_avogadro = 6.0221408e23;
const double mass_he = 0.004002602/n_avogadro; //kg
const double radius_he = 1.4e-10; //m
const double density_he = 0.1785; //kg/m3

//space between atoms surfaces and their square limit
const double r_ = pow(mass_he/(8*density_he), 1.0/3.0) - radius_he;


class simulate_n_particles{
  public:
    int N_, Nx;
    double sphere_radius;
    double x0, y0, z0;

    double* positions;
    double* velocities;

    ~simulate_n_particles(){
      //destructor
      delete[] positions;
      delete[] velocities;
    }

    simulate_n_particles(int N_, int Nx): N_(N_), Nx(Nx){
      //constructor
      sphere_radius = N_*(radius_he + r_)/2.0;
      x0 = radius_he + r_ - 4.0*sphere_radius, y0 = radius_he + r_ - 2.0*sphere_radius, z0 = radius_he + r_ - 2.0*sphere_radius;

      positions = new double[Nx*N_*N_*3];
      velocities = new double[Nx*N_*N_*3];

      cout << "N_ = "<< N_<<"\nR = " << sphere_radius << "\n";
    };


    int idx_grid(int i, int j, int k, int dim){
      return  (N_*(N_*i + j) + k)*3 + dim;
    }


    int set_initial_conditions(){
      int idx_x;
      for(int i = 0; i < Nx; i += 1){
        for(int j = 0; j < N_; j += 1){
          for(int k = 0; k < N_; k += 1){
            idx_x = idx_grid(i, j, k, 0);

            positions[idx_x] = x0 - 2*(radius_he + r_)*i;
            positions[idx_x + 1] = y0 + 2*(radius_he + r_)*j;
            positions[idx_x + 2] = z0 + 2*(radius_he + r_)*k;

            velocities[idx_x] = initial_speed_u;
            velocities[idx_x + 1] = 0.0;
            velocities[idx_x + 2] = 0.0;
//            save_position_file();
          }
        }
      }

      return 0;
    }


    int main(){
      set_initial_conditions();

      return 0;
    }


};
