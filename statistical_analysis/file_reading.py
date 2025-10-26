import os
import pandas as pd

def parse_meta(meta_str):
    data = meta_str.split(";")
    data = [x.split(" = ") for x in data]

    data_dict = {}
    for k, v in data:
        if str.strip(k) == "format_step_files":
            data_dict["format_step_files"] = v.split(", ")
            data_dict["steps_data"] = {x:[] for x in v.split(", ")}

            continue
        data_dict[str.strip(k)] = float(v)
    return data_dict


def data_Nx_N_(N_ : int, Nx : int):
    path_to_folder = os.path.join("..", "simulation_particles", "data", f"N_-{N_}_Nx-{Nx}")
    data_files = os.listdir(path_to_folder)
    data_files.remove("meta.txt")

    with open(os.path.join(path_to_folder, "meta.txt"), "r") as meta_data_file:
        meta_data_str = meta_data_file.read()
        data_dict = parse_meta(meta_data_str)

    n_steps = len(data_files)

    format_step_files = data_dict["format_step_files"]

    for n in range(n_steps):
        step_file_path = os.path.join(path_to_folder, f"step_{n}")

        with open(step_file_path, "r") as step_file:
            content = step_file.read()
            content = content.split("\n")[:-1]

        for particle_content in content:
            items = [float(x) for x in particle_content.split(", ")]

            for i in range(len(items)):
                variable_name = format_step_files[i]
                data_dict["steps_data"][variable_name].append(items[i])

    data_dict["steps_data"] = pd.DataFrame(data_dict["steps_data"])

    return data_dict


if __name__ == "__main__":
    data_Nx_N_(2, 1)
