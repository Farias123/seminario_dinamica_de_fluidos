from vpython import *

from file_reading import data_Nx_N_

def animate_Nx_N_(data_dict: dict, particles_to_animate : list) -> None:
  steps_data = data_dict["steps_data"]
  filtered_steps_data = steps_data[steps_data["idx"].isin(particles_to_animate)]

  R = data_dict["R"]
  scene.background = color.white
  central_sphere = sphere(pos=vector(0,0,0), radius=R, color=color.blue)
  scene.camera.pos = vector(5*R, 5*R, 5*R)
  scene.camera.axis = vector(-5*R, -5*R, -5*R)

  vpython_bodies = {}
  initial_positions = filtered_steps_data[filtered_steps_data["step"] == 0]
  for i in particles_to_animate:
    position_i = initial_positions.loc[initial_positions["idx"] == i, ["x", "y", "z"]]
    x, y, z = position_i.values[0]
    vpython_bodies[i] = sphere(pos=vector(x, y, z), make_trail=True,
                                 radius=data_dict['r'], retain=60, color=color.red)

  def update_positions(step: int) -> None:
    positions = filtered_steps_data[filtered_steps_data["step"] == step]

    for i in particles_to_animate:
      position_i = positions.loc[positions["idx"] == i, ["x", "y", "z"]]
      x, y, z = position_i.values[0]
      vpython_bodies[i].pos = vector(x, y, z)

  step_range = set(filtered_steps_data["step"])
  for step in step_range:
    rate(30)
    update_positions(step)


if __name__ == '__main__':
  Nx = 10
  N_ = 10

  data_dict = data_Nx_N_(Nx, N_)
  particles_to_animate = list(range(100))

  animate_Nx_N_(data_dict, particles_to_animate)