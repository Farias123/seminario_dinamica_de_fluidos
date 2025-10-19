#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

const double dt = 1e-15;
const double t_max = 1e-11; //s

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

    void save_position_file(int idx_x, double t){
      cout << "Baseado no número de partículas simuladas, se for menos que 10 mil, salva os dados de todas. Acima disso"
      << "sorteia 10 mil números aleatórios entre 0 e Nx*N_*N_. Se o índice da partícula (i + Nx*j + Nx*N_*k) é um dos"
      << "sorteados, será salvo no arquivo do passo de tempo. Esquema: criar uma pasta com nome tendo os N_ e Nx, para cada"
      << "passo criar arquivo com posições de todas as partículas";
    }

    void euler_update(int idx_x){
      positions[idx_x] += velocities[idx_x]*dt;
      positions[idx_x + 1] += velocities[idx_x + 1]*dt;
      positions[idx_x + 2] += velocities[idx_x + 2]*dt;

    //  velocities[idx_x]; Atualizar velocidade com colisao, para isso identificar colisao e depois resolver colisao

    }

    void update_and_save(double t){
      int idx_x; //x, y = idx_x + 1, z = idx_x + 2

      for(int i = 0; i < Nx; i += 1){
        for(int j = 0; j < N_; j += 1){
          for(int k = 0; k < N_; k += 1){
            idx_x = idx_grid(i, j, k, 0);

            euler_update(idx_x);
//            save_position_file(idx_x, t);

            cout << "N_ = " << N_ << ", idx x = " << idx_x << ", t = " << t << ", x = " << positions[idx_x] << ", y = " << positions[idx_x + 1] << ", z = " << positions[idx_x + 2] << "\n";
          }
        }
      }
    }

    void set_initial_conditions(){
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
    }

    void main(){
      set_initial_conditions();

      for(double t = 0; t < t_max; t += dt){
        update_and_save(t);
      }
    }


};
