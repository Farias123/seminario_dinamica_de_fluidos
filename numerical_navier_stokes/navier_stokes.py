import fenics as fe
from tqdm import tqdm  # optional
import matplotlib.pyplot as plt
# from mshr import *
import mshr
import numpy as np
from dolfin import *
import pyvista as pv


def navier_stoker_solver(temperature, sphere_radius, kinematic_viscosity, inlet_velocity, create_gif=False):
    L = 10 * sphere_radius  # Comprimento do túnel em X
    H = 4 * sphere_radius  # Comprimento do túnel em Y e Z
    U = inlet_velocity

    sphere_center = (L / 2, H / 2, H / 2)
    model_rank = 0
    sx, sy, sz = sphere_center

    N_POINTS_P_AXIS = 41
    KINEMATIC_VISCOSITY = kinematic_viscosity

    element_length = H / (N_POINTS_P_AXIS - 1)

    maximum_possible_time_step_length = (
            0.5 * element_length ** 2 / KINEMATIC_VISCOSITY
    )

    reynolds_number = U * H / KINEMATIC_VISCOSITY

    print(f"Reynolds Number: {reynolds_number}")
    print(maximum_possible_time_step_length)

    tol = 1
    t = 0.0
    T = 10.0  # Final time
    TIME_STEP_LENGTH = maximum_possible_time_step_length / 100  # Time step size
    N_TIME_STEPS = int(T / TIME_STEP_LENGTH)

    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], 0.0)

    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], L)

    class Bottom(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], 0.0)

    class Top(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[1], H)

    class caralho_z(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[2], 0.0)

    class piru_z(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[2], H)

    class SphereBoundary(SubDomain):
        def inside(self, x, on_boundary):
            dist = np.sqrt((x[0] - sx) ** 2 + (x[1] - sy) ** 2 + (x[2] - sz) ** 2)
            return on_boundary and near(dist, sphere_radius, tol)

    left = Left()
    top = Top()
    right = Right()
    bottom = Bottom()
    caralho_z = caralho_z()
    piru_z = piru_z()
    obstacle = SphereBoundary()

    box = mshr.Box(Point(0.0, 0.0, 0.0), Point(L, H, H))
    sphere = mshr.Sphere(Point(sx, sy, sz), sphere_radius)

    mesh = mshr.generate_mesh(box - sphere, N_POINTS_P_AXIS)

    domains = MeshFunction('size_t', mesh, mesh.topology().dim())

    domains.set_all(0)

    obstacle.mark(domains, 1)

    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    boundaries.set_all(0)

    obstacle.mark(boundaries, 1)

    left.mark(boundaries, 2)
    top.mark(boundaries, 3)
    right.mark(boundaries, 4)
    bottom.mark(boundaries, 5)
    caralho_z.mark(boundaries, 6)
    piru_z.mark(boundaries, 7)

    velocity_function_space = fe.VectorFunctionSpace(mesh, "Lagrange", 3)
    pressure_function_space = fe.FunctionSpace(mesh, "Lagrange", 1)

    u_trial = fe.TrialFunction(velocity_function_space)
    p_trial = fe.TrialFunction(pressure_function_space)
    v_test = fe.TestFunction(velocity_function_space)
    q_test = fe.TestFunction(pressure_function_space)

    object_boundary_condition = fe.DirichletBC(
        velocity_function_space,
        (0.0, 0.0, 0.0),
        boundaries,
        1
    )

    left_boundary_condition = fe.DirichletBC(
        velocity_function_space,
        (U, 0.0, 0.0),
        boundaries,
        2
    )

    top_boundary_condition = fe.DirichletBC(
        velocity_function_space,
        (0.0, 0.0, 0.0),
        boundaries,
        3
    )

    # right_boundary_condition = fe.DirichletBC(
    #     velocity_function_space,
    #     (0.0, 0.0, 0.0),
    #     boundaries,
    #     4
    # )

    bottom_boundary_condition = fe.DirichletBC(
        velocity_function_space,
        (0.0, 0.0, 0.0),
        boundaries,
        5
    )

    caralho_z_boundary_condition = fe.DirichletBC(
        velocity_function_space,
        (0.0, 0.0, 0.0),
        boundaries,
        6
    )

    piru_z_boundary_condition = fe.DirichletBC(
        velocity_function_space,
        (0.0, 0.0, 0.0),
        boundaries,
        7
    )

    pressure_left_boundary_condition = fe.DirichletBC(
        pressure_function_space,
        0.0,
        boundaries,
        2
    )

    pressure_right_boundary_condition = fe.DirichletBC(
        pressure_function_space,
        0.0,
        boundaries,
        4
    )

    velocity_boundary_conditions = [object_boundary_condition, left_boundary_condition, top_boundary_condition,
                                    bottom_boundary_condition, caralho_z_boundary_condition, piru_z_boundary_condition]
    pressure_boundary_conditions = [pressure_right_boundary_condition, pressure_left_boundary_condition]

    # Define the solution fields involved
    u_prev = fe.Function(velocity_function_space)
    u_tent = fe.Function(velocity_function_space)
    u_next = fe.Function(velocity_function_space)
    p_next = fe.Function(pressure_function_space)

    # Weak form of the momentum equation
    momentum_weak_form_residuum = (
            1.0 / TIME_STEP_LENGTH * fe.inner(u_trial - u_prev, v_test) * fe.dx
            +
            fe.inner(fe.grad(u_prev) * u_prev, v_test) * fe.dx
            +
            KINEMATIC_VISCOSITY * fe.inner(fe.grad(u_trial), fe.grad(v_test)) * fe.dx
    )
    momentum_weak_form_lhs = fe.lhs(momentum_weak_form_residuum)
    momentum_weak_form_rhs = fe.rhs(momentum_weak_form_residuum)

    # Weak form of the pressure poisson problem
    pressure_poisson_weak_form_lhs = fe.inner(fe.grad(p_trial), fe.grad(q_test)) * fe.dx
    pressure_poisson_weak_form_rhs = - 1.0 / TIME_STEP_LENGTH * fe.div(u_tent) * q_test * fe.dx

    # Weak form of the velocity update equation
    velocity_update_weak_form_lhs = fe.inner(u_trial, v_test) * fe.dx
    velocity_update_weak_form_rhs = (
            fe.inner(u_tent, v_test) * fe.dx
            -
            TIME_STEP_LENGTH * fe.inner(fe.grad(p_next), v_test) * fe.dx
    )

    # Pre-Compute the system matrices (because they do not greatly change)
    momentum_assembled_system_matrix = fe.assemble(momentum_weak_form_lhs)
    pressure_poisson_assembled_system_matrix = fe.assemble(pressure_poisson_weak_form_lhs)
    velocity_update_assembled_system_matrix = fe.assemble(velocity_update_weak_form_lhs)

    for t in tqdm(range(N_TIME_STEPS)):
        # (1) Solve for tentative velocity
        momentum_assembled_rhs = fe.assemble(momentum_weak_form_rhs)
        [bc.apply(momentum_assembled_system_matrix, momentum_assembled_rhs) for bc in velocity_boundary_conditions]
        fe.solve(
            momentum_assembled_system_matrix,
            u_tent.vector(),
            momentum_assembled_rhs,
            "gmres",
            "ilu",
        )

        # (2) Solve for the pressure
        pressure_poisson_assembled_rhs = fe.assemble(pressure_poisson_weak_form_rhs)
        [bc.apply(pressure_poisson_assembled_system_matrix, pressure_poisson_assembled_rhs) for bc in
         pressure_boundary_conditions]
        fe.solve(
            pressure_poisson_assembled_system_matrix,
            p_next.vector(),
            pressure_poisson_assembled_rhs,
            "gmres",
            "amg",
        )

        # (3) Correct the velocities to be incompressible
        [bc.apply(momentum_assembled_system_matrix, momentum_assembled_rhs) for bc in velocity_boundary_conditions]
        velocity_update_assembled_rhs = fe.assemble(velocity_update_weak_form_rhs)
        fe.solve(
            velocity_update_assembled_system_matrix,
            u_next.vector(),
            velocity_update_assembled_rhs,
            "gmres",
            "ilu",
        )

        # Advance in time
        u_prev.assign(u_next)

    # fe.plot(u_next, mode = "glyphs")
    # plt.show()
    u_file = fe.File("results/u_next.pvd")
    u_file << u_next

    return print("Pica boazona ai embaixo meu mano")


def generate_streamlines_data(velocity_archive, x_coord, y_coord, z_coord, sphere_radius, reescaling_factor):
    reader = pv.PVDReader(velocity_archive)

    mesh = reader.read()[0]

    vector_name = "f_38"

    # Criar sementes (exemplo: esfera de sementes)
    seeds = pv.Sphere(radius=sphere_radius, center=(5 * sphere_radius, 2 * sphere_radius, 2 * sphere_radius))

    # Gerar streamlines
    streamlines = mesh.streamlines_from_source(
        seeds,
        vectors=vector_name,
        max_steps=1000000,
        initial_step_length=0.0001,
        terminal_speed=1e-6,
        integrator_type=45  # Runge-Kutta 4/5
    )

    lines = streamlines.split_bodies()

    list_streamlines = []

    search_coordinates = np.array([x_coord, y_coord, z_coord])
    best_coordinates = np.array([0.0, 0.0, 0.0])
    searched_streamline = 0

    for sl in lines:
        pts = np.asarray(sl.points) * reescaling_factor
        list_streamlines.append(pts)

        for positions in range(len(pts)):
            position_vector = pts[positions]
            diference_vector = position_vector - search_coordinates
            diference_vector_modulus = np.sqrt(
                diference_vector[0] ** 2 + diference_vector[1] ** 2 + diference_vector[2] ** 2)
            algum_nome = position_vector - best_coordinates
            algum_nome_modulus = np.sqrt(algum_nome[0] ** 2 + algum_nome[1] ** 2 + algum_nome[2] ** 2)

            if diference_vector_modulus < algum_nome_modulus:
                best_coordinates = position_vector
                searched_streamline = sl

    return list_streamlines[int(searched_streamline)]


def create_streamline_graph(velocity_archive, sphere_radius):
    reader = pv.PVDReader(velocity_archive)

    mesh = reader.read()[0]

    vector_name = "f_38"

    # Criar sementes (exemplo: esfera de sementes)
    seeds = pv.Sphere(radius=4 * sphere_radius, center=(5 * sphere_radius, 2 * sphere_radius, 2 * sphere_radius))

    # Gerar streamlines
    streamlines = mesh.streamlines_from_source(
        seeds,
        vectors=vector_name,
        max_steps=10000,
        initial_step_length=0.01,
        terminal_speed=1e-6,
        integrator_type=45  # Runge-Kutta 4/5
    )

    # Plotar
    plotter = pv.Plotter()
    plotter.add_mesh(mesh, opacity=0.2)
    plotter.add_mesh(streamlines, line_width=2, color="blue")
    plotter.show()

    return print("Pica Top Bro")


if __name__ == "__main__":
    # navier_stoker_solver(temperature=273.5, sphere_radius=8.35, kinematic_viscosity=1, inlet_velocity=8)
    create_streamline_graph("results/u_next.pvd", sphere_radius=8.35)
    generate_streamlines_data(velocity_archive="results/u_next.pvd", x_coord=0.5, y_coord=0.5, z_coord=0.5,
                              sphere_radius=8.35, reescaling_factor=1e-7)