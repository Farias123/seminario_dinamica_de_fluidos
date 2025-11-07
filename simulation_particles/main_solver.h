#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <sstream>
#include <algorithm>
#include <unordered_set>
#include <filesystem>

using namespace std;

const double initial_speed_u = 1304.69; // m/s (T = 273.15 K)

const double dt = 1e-13; //s


//atom variables
const double n_avogadro = 6.0221408e23;
const double mass_he = 0.004002602/n_avogadro; //kg
const double radius_he = 1.4e-10; //m
const double density_he = 0.1785; //kg/m3

//space between atoms surfaces and their square limit
const double r_ = pow(mass_he/(8*density_he), 1.0/3.0) - radius_he;

vector<int> generate_random_particle_ids(int N_particles){
  srand(time(0));
  vector<int> selected_particles;

  if (N_particles <= 1000) {
      for (int i = 0; i < N_particles; i += 1) {
          selected_particles.push_back(i);
      }
  } else {
      unordered_set<int> unique_particles;

      while (unique_particles.size() < 1000) {
        int random_number = rand() % N_particles;
        unique_particles.insert(random_number);
      }
      selected_particles.assign(unique_particles.begin(), unique_particles.end());
  }

  return selected_particles;
}

class simulate_n_particles{
  public:
    int N_, Nx;
    double sphere_radius;
    double x0, y0, z0;
    const double t_max = 2*abs(x0 - 2*(radius_he + r_)*Nx)/initial_speed_u; //s

    double* positions;
    double* velocities;
    vector<int> saved_particles;
    string folder_name;

    ~simulate_n_particles(){
      //destructor
      delete[] positions;
      delete[] velocities;
    }

    simulate_n_particles(int N_, int Nx): N_(N_), Nx(Nx){
      //constructor
      sphere_radius = N_*(radius_he + r_)/2.0;
      x0 = radius_he + r_ - 4.0*sphere_radius, y0 = radius_he + r_ - 2.0*sphere_radius, z0 = radius_he + r_ - 2.0*sphere_radius;

      saved_particles = generate_random_particle_ids(Nx*N_*N_);

      stringstream temp_name;
      temp_name << "./data/N_-" << N_ << "_Nx-" << Nx;
      folder_name = temp_name.str();

      if (filesystem::exists(folder_name)) {
          filesystem::remove_all(folder_name);
      }

      create_save_folder();

      positions = new double[Nx*N_*N_*3];
      velocities = new double[Nx*N_*N_*3];

      cout << "N_ = "<< N_<<"\nR = " << sphere_radius << "\n";
    };

    void create_save_folder(){
      filesystem::create_directories(folder_name);

      stringstream file_name;
      file_name << folder_name << "/meta.txt";

      ofstream meta_file(file_name.str());

      meta_file << "r = " << radius_he << "; r_ = " << r_ << "; R = " << sphere_radius << "; dt = " << dt
      << "; U = "<< initial_speed_u << "; format_step_files = idx, t, x, y, z";
      meta_file.close();
    }

    void save_position_file(int idx_x, double t){
      if(find(saved_particles.begin(), saved_particles.end(), idx_x/3) != saved_particles.end()){ //particle is one of the random chosen to be saved
        stringstream file_name;
        file_name << folder_name << "/step_" << t/dt;

        ofstream step_file(file_name.str(), std::ios::app);
        step_file << idx_x/3 << ", " << t << ", " << positions[idx_x] << ", " << positions[idx_x + 1] << ", " << positions[idx_x + 2] << "\n";
        step_file.close();
      }
    }

    int idx_grid(int i, int j, int k, int dim){
      return  (N_*(N_*i + j) + k)*3 + dim;
    }

    void manage_collisions(int idx_x1){
      double x1, y1, z1, x2, y2, z2;
      double vx1, vy1, vz1, vx2, vy2, vz2;
      double alpha, beta, d, intersection, nx, ny, nz;
      int idx_x2;

      x1 = positions[idx_x1], y1 = positions[idx_x1 + 1], z1 = positions[idx_x1 + 2];
      vx1 = velocities[idx_x1], vy1 = velocities[idx_x1 + 1], vz1 = velocities[idx_x1 + 2];
      d = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2)); //distance to origin (sphere with radius = R)

      //collision with central sphere
      if(d < radius_he + sphere_radius){
        nx = - x1/d, ny = -y1/d, nz = - z1/d;

        intersection = radius_he + sphere_radius - d;

        // position correction before changing velocity
        positions[idx_x1] = x1 - intersection*nx;
        positions[idx_x1 + 1] = y1 - intersection*ny;
        positions[idx_x1 + 2] = z1 - intersection*nz;

        velocities[idx_x1] = vx1 - 2*nx*(vx1*nx + vy1*ny + vz1*nz);
        velocities[idx_x1 + 1] = vy1 - 2*ny*(vx1*nx + vy1*ny + vz1*nz);
        velocities[idx_x1 + 2] = vz1 - 2*nz*(vx1*nx + vy1*ny + vz1*nz);
      }

      // collisions between particles
      for(int i = 0; i < Nx; i += 1){
        for(int j = 0; j < N_; j += 1){
          for(int k = 0; k < N_; k += 1){
            idx_x2 = idx_grid(i, j, k, 0);

            if(idx_x1 != idx_x2){
              x1 = positions[idx_x1], y1 = positions[idx_x1 + 1], z1 = positions[idx_x1 + 2];
              x2 = positions[idx_x2], y2 = positions[idx_x2 + 1], z2 = positions[idx_x2 + 2];
              d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2)); //distance between particles

              if(d < 2*radius_he){
                nx = (x2 - x1)/d, ny = (y2 - y1)/d, nz = (z2 - z1)/d;

                intersection = 2*radius_he - d;

                // position correction before changing velocity
                positions[idx_x1] = x1 - intersection/2*nx;
                positions[idx_x1 + 1] = y1 - intersection/2*ny;
                positions[idx_x1 + 2] = z1 - intersection/2*nz;

                positions[idx_x2] = x2 + intersection/2*nx;
                positions[idx_x2 + 1] = y2 + intersection/2*ny;
                positions[idx_x2 + 2] = z2 + intersection/2*nz;

                // velocities update
                vx1 = velocities[idx_x1], vy1 = velocities[idx_x1 + 1], vz1 = velocities[idx_x1 + 2];
                vx2 = velocities[idx_x2], vy2 = velocities[idx_x2 + 1], vz2 = velocities[idx_x2 + 2];

                alpha = (vx2 - vx1)*nx + (vy2 - vy1)*ny + (vz2 - vz1)*nz; // factors come from vector projection
                beta = (vx1 - vx2)*nx + (vy1 - vy2)*ny + (vz1 - vz2)*nz;

                velocities[idx_x1] = vx1 + nx*alpha;
                velocities[idx_x1 + 1] = vy1 + ny*alpha;
                velocities[idx_x1 + 2] = vz1 + nz*alpha;

                velocities[idx_x2] = vx2 + nx*beta;
                velocities[idx_x2 + 1] = vy2 + ny*beta;
                velocities[idx_x2 + 2] = vz2 + nz*beta;
              }

            }

          }
        }
      }

    }

    void euler_update(int idx_x){
      positions[idx_x] += velocities[idx_x]*dt;
      positions[idx_x + 1] += velocities[idx_x + 1]*dt;
      positions[idx_x + 2] += velocities[idx_x + 2]*dt;

      manage_collisions(idx_x);
    }

    void update_and_save(double t){
      int idx_x; //x, y = idx_x + 1, z = idx_x + 2

      for(int i = 0; i < Nx; i += 1){
        for(int j = 0; j < N_; j += 1){
          for(int k = 0; k < N_; k += 1){
            idx_x = idx_grid(i, j, k, 0);

            euler_update(idx_x);
            save_position_file(idx_x, t);
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
            save_position_file(idx_x, 0.0);
          }
        }
      }
    }

    void main(){
      set_initial_conditions();

      for(double t = dt; t < t_max; t += dt){
        update_and_save(t);
      }
//      todo: dividir espaco em caixas para otimizar deteccao de colisao; calcular tempo aproximado para particula atravessar esfera, dt e R.

    }
};
