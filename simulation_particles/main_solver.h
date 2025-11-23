#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <filesystem>
#include <chrono>

using namespace std;
using namespace std::chrono;

const double initial_speed_u = 1304.69; // m/s (T = 273.15 K)

const double dt = 1e-13; //s
//const double sphere_radius = 8.34816651e-7; //m

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
    int N_, Nx, N;
    double sphere_radius;
    double x0, y0, z0;
    const double t_max = 2*abs(x0 - 2*(radius_he + r_)*Nx)/initial_speed_u; //s

    double* x;
    double* y;
    double* z;
    double* vx;
    double* vy;
    double* vz;

    vector<int> saved_particles;
    string folder_name;

    high_resolution_clock::time_point simulation_start_time = high_resolution_clock::now();

    ~simulate_n_particles(){
      //destructor
      delete[] x;
      delete[] y;
      delete[] z;
      delete[] vx;
      delete[] vy;
      delete[] vz;
    }

    simulate_n_particles(int N_, int Nx): N_(N_), Nx(Nx), N(Nx*N_*N_){
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

      filesystem::create_directories(folder_name);
      x = new double[Nx*N_*N_];
      y = new double[Nx*N_*N_];
      z = new double[Nx*N_*N_];
      vx = new double[Nx*N_*N_];
      vy = new double[Nx*N_*N_];
      vz = new double[Nx*N_*N_];

      cout << "N_ = "<< N_<<"\nR = " << sphere_radius << "\n";
    };

    void create_meta_file(){
      stringstream file_name;

      file_name << folder_name << "/meta.txt";

      ofstream meta_file(file_name.str());

      auto simulation_end_time = high_resolution_clock::now();
      auto simulation_duration = duration_cast<milliseconds>(simulation_end_time - simulation_start_time);

      meta_file << "r = " << radius_he << "; r_ = " << r_ << "; R = " << sphere_radius << "; dt = " << dt
      << "; U = "<< initial_speed_u <<"; simulation_time(ms) = "<< to_string(simulation_duration.count()) << "; format_step_files = idx, x, y, z";
      meta_file.close();
    }

    void save_position_file(double t){
      stringstream file_name;
      file_name << folder_name << "/step_" << t/dt;
      ofstream step_file(file_name.str(), std::ios::app);

      for(int idx_particle = 0; idx_particle < N; idx_particle += 1){
        if(find(saved_particles.begin(), saved_particles.end(), idx_particle) != saved_particles.end()){ //particle is one of the random chosen to be saved
          step_file << idx_particle << ", " << x[idx_particle] << ", " << y[idx_particle] << ", " << z[idx_particle] << "\n";
        }
      }
      step_file.close();
    }

    int idx_grid(int i, int j, int k){
      return  N_*(N_*i + j) + k;
    }

    void manage_collisions(int idx_p1){
      double x1, y1, z1, x2, y2, z2;
      double vx1, vy1, vz1, vx2, vy2, vz2;
      double alpha, beta, d, intersection, nx, ny, nz;

      x1 = x[idx_p1], y1 = y[idx_p1], z1 = z[idx_p1];
      vx1 = vx[idx_p1], vy1 = vy[idx_p1], vz1 = vz[idx_p1];

      d = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2)); //distance to origin (sphere with radius = R)

      //collision with central sphere
      if(d < radius_he + sphere_radius){
        nx = - x1/d, ny = -y1/d, nz = - z1/d;

        intersection = radius_he + sphere_radius - d;

        // position correction before changing velocity
        x[idx_p1] = x1 - intersection*nx;
        y[idx_p1] = y1 - intersection*ny;
        z[idx_p1] = z1 - intersection*nz;

        vx[idx_p1] = vx1 - 2*nx*(vx1*nx + vy1*ny + vz1*nz);
        vy[idx_p1] = vy1 - 2*ny*(vx1*nx + vy1*ny + vz1*nz);
        vz[idx_p1] = vz1 - 2*nz*(vx1*nx + vy1*ny + vz1*nz);
      }

      // collisions between particles
      for(int idx_p2 = 0; idx_p2 < Nx; idx_p2 += 1){
        if(idx_p1 != idx_p2){
          x2 = x[idx_p2], y2 = y[idx_p2], z2 = z[idx_p2];
          d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2)); //distance between particles

          if(d < 2*radius_he){
            nx = (x2 - x1)/d, ny = (y2 - y1)/d, nz = (z2 - z1)/d;

            intersection = 2*radius_he - d;

            // position correction before changing velocity
            x[idx_p1] = x1 - intersection/2*nx;
            y[idx_p1] = y1 - intersection/2*ny;
            z[idx_p1] = z1 - intersection/2*nz;

            x[idx_p2] = x2 + intersection/2*nx;
            y[idx_p2] = y2 + intersection/2*ny;
            z[idx_p2] = z2 + intersection/2*nz;

            // velocities update
            vx2 = vx[idx_p2], vy2 = vy[idx_p2], vz2 = vz[idx_p2];

            alpha = (vx2 - vx1)*nx + (vy2 - vy1)*ny + (vz2 - vz1)*nz; // factors come from vector projection
            beta = (vx1 - vx2)*nx + (vy1 - vy2)*ny + (vz1 - vz2)*nz;

            vx[idx_p1] = vx1 + nx*alpha;
            vy[idx_p1] = vy1 + ny*alpha;
            vz[idx_p1] = vz1 + nz*alpha;

            vx[idx_p2] = vx2 + nx*beta;
            vy[idx_p2] = vy2 + ny*beta;
            vz[idx_p2] = vz2 + nz*beta;
          }
        }
      }

    }

    void euler_update(int idx_particle){
      x[idx_particle] += vx[idx_particle]*dt;
      y[idx_particle] += vy[idx_particle]*dt;
      z[idx_particle] += vz[idx_particle]*dt;

      manage_collisions(idx_particle);
    }

    void update_and_save(double t){
      for(int idx_particle = 0; idx_particle < N; idx_particle += 1){
        euler_update(idx_particle);
      }
      save_position_file(t);
    }

    void set_initial_conditions(){
      int idx_particle;
      for(int i = 0; i < Nx; i += 1){
        for(int j = 0; j < N_; j += 1){
          for(int k = 0; k < N_; k += 1){
            idx_particle = idx_grid(i, j, k);

            x[idx_particle] = x0 - 2*(radius_he + r_)*i;
            y[idx_particle] = y0 + 2*(radius_he + r_)*j;
            z[idx_particle] = z0 + 2*(radius_he + r_)*k;

            vx[idx_particle] = initial_speed_u;
            vy[idx_particle] = 0.0;
            vz[idx_particle] = 0.0;
          }
        }
      }
      save_position_file(0.0);
    }

    void main(){
      set_initial_conditions();

      for(double t = dt; t < t_max; t += dt){
        update_and_save(t);
      }

      create_meta_file();
//      todo: dividir espaco em caixas para otimizar deteccao de colisao;

    }
};
