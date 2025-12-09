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

const double n_avogadro = 6.0221408e23;
const double boltzmann_constant = 1.380649e-23; //m2* kg/(K * s2)

//atom variables - water molecule
const double mass_atom = 0.0180153/n_avogadro; //kg
const double radius_atom = 2.8e-10; //m source: https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=6&id=103723

//boundary conditions
const double temperature = 273.15; //K
const double initial_speed_u = sqrt(3*boltzmann_constant*temperature / mass_atom); //m/s (T = 273.15 K)
const double asymptotic_U = 8.35e-7; //m/s

const double dt = 1e-11; //s
const double time_to_equilibrium = 5000.0*dt; //s
const double sphere_radius = 8.35e-7; //m
const double x_sphere = 5.0*sphere_radius;
const double y_sphere = 2.0*sphere_radius;
const double z_sphere = 2.0*sphere_radius;

const double sigma_v = sqrt(boltzmann_constant*temperature / mass_atom); //sqrt(kb*T/m)

vector<int> generate_random_particle_ids(int N_particles){
  srand(time(0));
  vector<int> selected_particles;
  int n_particles_to_save = 10000;

  if (N_particles <= n_particles_to_save) {
      for (int i = 0; i < N_particles; i += 1) {
          selected_particles.push_back(i);
      }
  } else {
      unordered_set<int> unique_particles;

      while (unique_particles.size() < n_particles_to_save) {
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
    bool save_positions;
    int number_of_collisions_step;

    const double x_a, x_b, y_a, y_b, z_a, z_b;
    const double t_max; //s
    const double Lcell;

    const int Cx, Cy, Cz;
    const long int N_cells_total;

    double* x;
    double* y;
    double* z;
    double* vx;
    double* vy;
    double* vz;

    // linked list
    int* head;
    int* linked_list;

    vector<int> saved_particles;
    string folder_name;
    string collisions_step_filename;

    high_resolution_clock::time_point simulation_start_time = high_resolution_clock::now();

    ~simulate_n_particles(){
      //destructor
      delete[] x;
      delete[] y;
      delete[] z;
      delete[] vx;
      delete[] vy;
      delete[] vz;

      delete[] head;
      delete[] linked_list;
    }

    simulate_n_particles(int N, bool save_positions): N(N), save_positions(save_positions),
    t_max(100*dt), Lcell(15000*radius_atom), x_a(0.0), x_b(4.0*sphere_radius), y_a(0.0),
    y_b(4.0*sphere_radius), z_a(0.0), z_b(4.0*sphere_radius), Cx(2*(x_sphere - x_a) / Lcell), Cy((y_b - y_a) / Lcell),
    Cz((z_b - z_a) / Lcell), N_cells_total(Cx*Cy*Cz){
      //constructor

      saved_particles = generate_random_particle_ids(N);

      stringstream temp_name;
      temp_name << "./data/N-" << N;
      folder_name = temp_name.str();

      temp_name.str("");
      temp_name.clear();
      temp_name << folder_name << "./data/" << N << "_number_collisions_step.csv";
      collisions_step_filename = temp_name.str();

      if (filesystem::exists(folder_name) && (save_positions == true)) {
          filesystem::remove_all(folder_name);
      }

      filesystem::create_directories(folder_name);
      filesystem::remove(collisions_step_filename);

      x = new double[N];
      y = new double[N];
      z = new double[N];
      vx = new double[N];
      vy = new double[N];
      vz = new double[N];

      head = new int[N_cells_total];
      linked_list = new int[N];

      cout << "N = "<< N << "\n";
    };

    void create_meta_file(){
      stringstream file_name;

      file_name << folder_name << "/meta.txt";

      ofstream meta_file(file_name.str());

      auto simulation_end_time = high_resolution_clock::now();
      auto simulation_duration = duration_cast<milliseconds>(simulation_end_time - simulation_start_time);

      meta_file << "r = " << radius_atom << "; R = " << sphere_radius << "; dt = " << dt
      << "; U = "<< asymptotic_U <<"; simulation_time(ms) = "<< to_string(simulation_duration.count()) << "; format_step_files = idx, x, y, z, vx, vy, vz";
      meta_file.close();
    }

    void save_position_file(double t){
      stringstream file_name;
      file_name << folder_name << "/step_" << t/dt;
      ofstream step_file(file_name.str(), std::ios::app);

      for(int idx_particle : saved_particles){ // save selected particles
        step_file << idx_particle << ", " << x[idx_particle] << ", " << y[idx_particle] << ", " << z[idx_particle] <<
        ", " << vx[idx_particle] << ", " << vy[idx_particle] << ", " << vz[idx_particle] << "\n";
      }
      step_file.close();
    }

    int idx_cell_grid(int i, int j, int k){
      return Cz*(Cy*i + j) + k;
    }

    void reset_particle_cells(){
      int idx_cell;

      for(idx_cell = 0; idx_cell < N_cells_total; idx_cell += 1){
        head[idx_cell] = -1;
      }
    }

    void particle_collisions_linked_list(){
      double x1, y1, z1, x2, y2, z2;
      double vx1, vy1, vz1, vx2, vy2, vz2;
      double alpha, beta, d, intersection, nx, ny, nz;
      int idx_p1, idx_p2;
      int cx,cy,cz,c;
      int ccx,ccy,ccz;
      int cxviz,cyviz,czviz,cviz;
      int qual[3];

      reset_particle_cells();

      for(int i = 0; i < N; i += 1){
        qual[0] = (x[i] - x_a)/Lcell;
        qual[1] = (y[i] - y_a)/Lcell;
        qual[2] = (z[i] - z_a)/Lcell;
        if((qual[0] < 0) || (qual[0] > Cx) || (qual[1] < 0) || (qual[1] > Cy) || (qual[2] < 0) || (qual[2] > Cz))
          continue;

        c = idx_cell_grid(qual[0], qual[1], qual[2]);

        linked_list[i] = head[c];
        //a última será a head / inicial
        head[c] = i;
      }

      for(cx = 1; cx < Cx - 1; cx += 1)
      for(cy = 1; cy < Cy - 1; cy += 1)
      for(cz = 1; cz < Cz - 1; cz += 1){
        c = idx_cell_grid(cx, cy, cz);

        for(ccx = - 1; ccx <= 1; ccx += 1)
        for(ccy = - 1; ccy <= 1; ccy += 1)
        for(ccz = - 1; ccz <= 1; ccz += 1){
          cxviz = cx + ccx;
          cyviz = cy + ccy;
          czviz = cz + ccz;

          cviz = idx_cell_grid(cxviz, cyviz, czviz);

          idx_p1 = head[c];

          while(idx_p1 >= 0){
            x1 = x[idx_p1], y1 = y[idx_p1], z1 = z[idx_p1];
            vx1 = vx[idx_p1], vy1 = vy[idx_p1], vz1 = vz[idx_p1];

            idx_p2 = head[cviz];

            while(idx_p2 >= 0){

              if(idx_p1 != idx_p2){
                x2 = x[idx_p2], y2 = y[idx_p2], z2 = z[idx_p2];

                d = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2)); //distance between particles

                if(d < 2*radius_atom){
                  nx = (x2 - x1)/d, ny = (y2 - y1)/d, nz = (z2 - z1)/d;
                  intersection = 2*radius_atom - d;

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

              idx_p2 = linked_list[idx_p2];
            }
            idx_p1 = linked_list[idx_p1];
          }
        }
      }
    }

    void sphere_and_wall_collisions(int idx_p1){
      double x1, y1, z1;
      double vx1, vy1, vz1;
      double d, intersection, nx, ny, nz;

      x1 = x[idx_p1], y1 = y[idx_p1], z1 = z[idx_p1];
      vx1 = vx[idx_p1], vy1 = vy[idx_p1], vz1 = vz[idx_p1];

      d = sqrt(pow(x1 - x_sphere, 2) + pow(y1 - y_sphere, 2) + pow(z1 - z_sphere, 2)); //distance to origin (sphere with radius = R)

      //collision with central sphere
      if(d < radius_atom + sphere_radius){
        nx = - (x1 - x_sphere)/d, ny = - (y1 - y_sphere)/d, nz = - (z1 - z_sphere)/d;

        intersection = radius_atom + sphere_radius - d;

        // position correction before changing velocity
        x[idx_p1] = x1 - intersection*nx;
        y[idx_p1] = y1 - intersection*ny;
        z[idx_p1] = z1 - intersection*nz;

        vx[idx_p1] = vx1 - 2*nx*(vx1*nx + vy1*ny + vz1*nz);
        vy[idx_p1] = vy1 - 2*ny*(vx1*nx + vy1*ny + vz1*nz);
        vz[idx_p1] = vz1 - 2*nz*(vx1*nx + vy1*ny + vz1*nz);

        number_of_collisions_step += 1;
      }

      //collision walls
      d = y_b - y[idx_p1]; //right
      if(d < radius_atom){
        y[idx_p1] -= (radius_atom - d);
        vy[idx_p1] = - vy[idx_p1];
      }

      d = - (y_a - y[idx_p1]); //left
      if(d < radius_atom){
        y[idx_p1] += (radius_atom - d);
        vy[idx_p1] = - vy[idx_p1];
      }

      d = z_b - z[idx_p1]; //up
      if(d < radius_atom){
        z[idx_p1] -= (radius_atom - d);
        vz[idx_p1] = - vz[idx_p1];
      }

      d = - (z_a - z[idx_p1]); //down
      if(d < radius_atom){
        z[idx_p1] += (radius_atom - d);
        vz[idx_p1] = - vz[idx_p1];
      }

      //check if it is out of the space
      if((x[idx_p1] < x_a) || (x[idx_p1] > x_a + 2*(x_sphere - x_a))){
        assign_random_valid_position(idx_p1, false);
        assign_random_velocity(idx_p1);
      }

    }

    void euler_update(){
      for(int idx_particle = 0; idx_particle < N; idx_particle += 1){
        x[idx_particle] += vx[idx_particle]*dt;
        y[idx_particle] += vy[idx_particle]*dt;
        z[idx_particle] += vz[idx_particle]*dt;

        sphere_and_wall_collisions(idx_particle);
      }
    }

    void update(){
      euler_update();
      particle_collisions_linked_list();
    }

    void save_number_of_sphere_collisions_to_file(double t){
      ofstream sphere_collisions_file(collisions_step_filename, std::ios::app);
      sphere_collisions_file << t << ", " << number_of_collisions_step << "\n";

      sphere_collisions_file.close();
    }

    void assign_random_valid_position(int idx_p1, bool is_first_generation){
      double temp_x, temp_y, temp_z;
      double x2, y2, z2;

      int idx_p2;
      int cx,cy,cz, c;
      int ccx,ccy,ccz;
      int cxviz,cyviz,czviz,cviz;
      int qual[3];

      uniform_real_distribution<double> x_distribution(x_a + radius_atom, x_b - radius_atom);
      uniform_real_distribution<double> y_distribution(y_a + radius_atom, y_b - radius_atom);
      uniform_real_distribution<double> z_distribution(z_a + radius_atom, z_b - radius_atom);

      double distance_to_existing_particle;

      static std::random_device rd; //numeros aleatorios de alta precisao
      static std::mt19937 generator(rd());

      bool pos_is_valid = false;

      reset_particle_cells();

      int particles_already_created = N;
      if(is_first_generation){
        particles_already_created = idx_p1;
      }

      for(int i = 0; i < particles_already_created; i += 1){
        qual[0] = (x[i] - x_a)/Lcell;
        qual[1] = (y[i] - y_a)/Lcell;
        qual[2] = (z[i] - z_a)/Lcell;
        if((qual[0] < 0) || (qual[0] >= Cx) || (qual[1] < 0) || (qual[1] >= Cy) || (qual[2] < 0) || (qual[2] >= Cz))
          continue;

        c = idx_cell_grid(qual[0], qual[1], qual[2]);

        linked_list[i] = head[c];
        //a última será a head / inicial
        head[c] = i;
      }

      while(pos_is_valid == false){
        temp_x = x_distribution(generator);
        temp_y = y_distribution(generator);
        temp_z = z_distribution(generator);

        cx = (temp_x - x_a)/Lcell;
        cy = (temp_y - y_a)/Lcell;
        cz = (temp_z - z_a)/Lcell;

        pos_is_valid = true;

        for(ccx = - 1; ccx <= 1 && pos_is_valid; ccx += 1)
        for(ccy = - 1; ccy <= 1 && pos_is_valid; ccy += 1)
        for(ccz = - 1; ccz <= 1 && pos_is_valid; ccz += 1){
          cxviz = cx + ccx;
          cyviz = cy + ccy;
          czviz = cz + ccz;

          if((cxviz < 0) || (cxviz >= Cx) || (cyviz < 0) || (cyviz >= Cy) || (czviz < 0) || (czviz >= Cz))
            continue;

          cviz = idx_cell_grid(cxviz, cyviz, czviz);

          idx_p2 = head[cviz];

          while(idx_p2 >= 0){
            x2 = x[idx_p2], y2 = y[idx_p2], z2 = z[idx_p2];

            distance_to_existing_particle = sqrt(pow(x2 - temp_x, 2) + pow(y2 - temp_y, 2) + pow(z2 - temp_z, 2));
            if(distance_to_existing_particle < 2*radius_atom){
              pos_is_valid = false;
              break;
            }

            idx_p2 = linked_list[idx_p2];
          }
        }
      }
      x[idx_p1] = temp_x;
      y[idx_p1] = temp_y;
      z[idx_p1] = temp_z;
    }

    void assign_random_velocity(int idx_p1){
      static std::random_device rd;
      static std::mt19937 generator(rd());
      double temp_vx = -1;

      normal_distribution maxwell_boltzmann_dist{0.0, sigma_v};

      //dist normal (vx soma U)
      while(temp_vx < 0){
        temp_vx = maxwell_boltzmann_dist(generator);
      }
      vx[idx_p1] = temp_vx + asymptotic_U;
      vy[idx_p1] = maxwell_boltzmann_dist(generator);
      vz[idx_p1] = maxwell_boltzmann_dist(generator);
    }

    void set_initial_conditions(){
      bool is_first_generation = true;

      for(int idx_p1 = 0; idx_p1 < N; idx_p1 += 1){
        assign_random_valid_position(idx_p1, is_first_generation);
        assign_random_velocity(idx_p1);
      }

      if(save_positions == true)
        save_position_file(0.0);
    }

    void main(){
      set_initial_conditions();
      for(double t = 0.0; t < time_to_equilibrium; t += dt)
        update();

      number_of_collisions_step = 0;
      save_number_of_sphere_collisions_to_file(0.0);

      for(double t = dt; t < t_max; t += dt){
        number_of_collisions_step = 0;
        update();
        save_number_of_sphere_collisions_to_file(t);

        if(save_positions == true)
          save_position_file(t);
      }

      create_meta_file();
    }
};
