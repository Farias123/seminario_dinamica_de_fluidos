import subprocess
import os
import shutil

from scipy.signal import step2

from statistical_analysis.file_reading import data_N
from numerical_navier_stokes.navier_stokes import extract_velocity_field

import pandas as pd
import numpy as np

def find_most_similar_row(row, velocity_field):
    x, y, z = row[["x", "y", "z"]]

    pos_fields = [t[0] for t in velocity_field]

    n = len(pos_fields)
    proximity = np.zeros(n)

    for i in range(n):
        proximity[i] = (pos_fields[i][0] - x)**2 + (pos_fields[i][1] - y)**2 + (pos_fields[i][2] - z)**2

    min_index = np.argmin(proximity)

    return velocity_field[min_index]


def full_process(n_end):
    subprocess.run(["g++", "../simulation_particles/main.cpp", "-o",  "./main"], capture_output=True, text=True)

    folder_path = os.path.join("..", "statistical_analysis", "data")
    content = os.listdir(folder_path)

    file_path = os.path.join(folder_path, "residue.csv")

    if os.path.exists(file_path):
        content = pd.read_csv(file_path, header=None)
        n_start = content.iloc[:,0].max()
    else:
        n_start = 100000

    velocity_archive = os.path.join("..", "numerical_navier_stokes", "results", "u_next.pvd")
    velocity_field = extract_velocity_field(velocity_archive, 1e-7)
    for n in range(n_start, n_end, 10000):
        params = ["./main", str(n)]
        output = subprocess.run(params, capture_output=True, text=True)

        path_to_folder = os.path.join("..", "statistical_analysis", "data", f"N-{n}")

        data_dict = data_N(path_to_folder)
        steps_data = data_dict["steps_data"]

        residue = 0
        for i in range(len(steps_data)):
            row = steps_data.iloc[i]

            most_similar_row = find_most_similar_row(row, velocity_field)

            row_velocities = row[["vx", "vy", "vz"]]
            residue += np.linalg.norm(most_similar_row[1] - row[["vx", "vy", "vz"]])

        if not os.path.exists(file_path):
            with open(file_path, "w") as file:
                file.write("")

        with open(file_path, "a") as file:
            residue_string = f"{n}, {residue}\n"
            file.write(residue_string)

        shutil.rmtree(path_to_folder)


if __name__ == "__main__":
    n_end = 1000000
    full_process(n_end)