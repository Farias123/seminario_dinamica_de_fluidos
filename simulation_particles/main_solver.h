#include <cmath>
#include <array>
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
const double sphere_radius = 8.34816651e-7; //m

//atom variables
const double n_avogadro = 6.0221408e23;
const double mass_he = 0.004002602/n_avogadro; //kg
const double radius_he = 1.4e-10; //m
const double density_he = 0.1785; //kg/m3

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
    int N;
    double x_a, x_b, y_a, y_b, z_a, z_b;
    const double t_max; //s

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

    simulate_n_particles(int N): N(N), t_max(12*sphere_radius/(initial_speed_u)){
      //constructor
      x_a = -6.0*sphere_radius, x_b = -2.0*sphere_radius;
      y_a = -2.0*sphere_radius, y_b = 2.0*sphere_radius;
      z_a = -2.0*sphere_radius, z_b = 2.0*sphere_radius;

      saved_particles = generate_random_particle_ids(N);

      stringstream temp_name;
      temp_name << "./data/N-" << N;
      folder_name = temp_name.str();

      if (filesystem::exists(folder_name)) {
          filesystem::remove_all(folder_name);
      }

      filesystem::create_directories(folder_name);
      x = new double[N];
      y = new double[N];
      z = new double[N];
      vx = new double[N];
      vy = new double[N];
      vz = new double[N];

      cout << "N = "<< N << "\n";
    };

    void create_meta_file(){
      stringstream file_name;

      file_name << folder_name << "/meta.txt";

      ofstream meta_file(file_name.str());

      auto simulation_end_time = high_resolution_clock::now();
      auto simulation_duration = duration_cast<milliseconds>(simulation_end_time - simulation_start_time);

      meta_file << "r = " << radius_he << "; R = " << sphere_radius << "; dt = " << dt
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
      for(int idx_p2 = 0; idx_p2 < N; idx_p2 += 1){
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

    array<double, 3> generate_random_pos(){
      double temp_x, temp_y, temp_z;

      default_random_engine generator(time(0));
      uniform_real_distribution<double> x_distribution(x_a, x_b);
      uniform_real_distribution<double> y_distribution(y_a, y_b);
      uniform_real_distribution<double> z_distribution(z_a, z_b);

      temp_x = x_distribution(generator);
      temp_y = y_distribution(generator);
      temp_z = z_distribution(generator);

      return {temp_x, temp_y, temp_z};
    }

    void set_initial_conditions(){
      double temp_x, temp_y, temp_z;
      double x2, y2, z2;

      double distance_to_existing_particle;
      array<double, 3> temp_pos;

      bool pos_is_valid;

      for(int idx_p1 = 0; idx_p1 < N; idx_p1 += 1){

        pos_is_valid = false;
        while(pos_is_valid == false){
          temp_pos = generate_random_pos();
          temp_x = temp_pos[0];
          temp_y = temp_pos[1];
          temp_z = temp_pos[2];

          pos_is_valid = true;

          for(int idx_p2 = 0; idx_p2 < idx_p1; idx_p2 += 1){
            x2 = x[idx_p2], y2 = y[idx_p2], z2 = z[idx_p2];

            distance_to_existing_particle = sqrt(pow(x2 - temp_x, 2) + pow(y2 - temp_y, 2) + pow(z2 - temp_z, 2));

            if(distance_to_existing_particle < 2*radius_he){
              pos_is_valid = false;
              break;
            }
          }

          if (pos_is_valid){
            x[idx_p1] = temp_x;
            y[idx_p1] = temp_y;
            z[idx_p1] = temp_z;
          }
        }

        vx[idx_p1] = initial_speed_u;
        vy[idx_p1] = 0.0;
        vz[idx_p1] = 0.0;
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
