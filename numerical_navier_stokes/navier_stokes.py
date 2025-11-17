import matplotlib.pyplot as plt
import numpy as np
import fipy as fp
from fipy import *
import pyvista as pv
from tqdm import tqdm


def main():

    #Boundary
    Lx = Ly = Lz = 2

    #Nodes
    nx = ny = nz = 41

    #Assintotic Velocity

    U = 10.0

    #Build the grid

    mesh = Grid3D(Lx = Lx, Ly = Ly, Lz = Lz, nx = nx, ny = ny, nz = nz)

    vx = CellVariable(name = "x-velocity", mesh=mesh, hasOld= True)
    vy = CellVariable(name="y-velocity", mesh=mesh, hasOld=True)
    vz = CellVariable(name="z-velocity", mesh=mesh, hasOld=True)

    v = FaceVariable(name = "velocity", mesh = mesh , rank = 1)

    p = CellVariable(name = "pressure", mesh = mesh, hasOld= True)

    #Boundary Conditions

    #top
    vy.constrain(0.0,where = mesh.facesTop)
    p.faceGrad.dot([0,1,0]).constrain(0.0, where=mesh.facesTop)

    #left
    vx.constrain(U, where = mesh.facesLeft)
    p.constrain(0.0, where=mesh.facesLeft)

    #right
    p.constrain(0.0, where=mesh.facesRight)

    #bottom
    vy.constrain(0.0, where = mesh.facesBottom)
    p.faceGrad.dot([0,-1,0]).constrain(0.0, where=mesh.facesBottom)

    #front
    vz.constrain(0.0, where = mesh.facesFront)
    p.faceGrad.dot([0,0,-1]).constrain(0.0, where=mesh.facesFront)

    #back
    vz.constrain(0.0, where=mesh.facesBack)
    p.faceGrad.dot([0,0,1]).constrain(0.0, where=mesh.facesBack)


    nu = 0.1 #kinematic_viscosity
    rho = 1.0 #Densidade
    dt = 0.01

    #Equations

    eqv_i = (TransientTerm() == DiffusionTerm(nu) - ConvectionTerm(v))
    eq_p = (DiffusionTerm(coeff=1.) == (rho/dt)*v.divergence)

    steps = 0

    for __ in tqdm(range(steps)):
        vx.updateOld()
        vy.updateOld()
        vz.updateOld()

        v[0, :] = vx.faceValue
        v[1, :] = vy.faceValue
        v[2, :] = vz.faceValue

        res0 = eqv_i.solve(var=vx, dt=dt)
        res1 = eqv_i.solve(var=vy, dt=dt)
        res2 = eqv_i.solve(var=vz, dt=dt)

        res_p = eq_p.solve(var=p, dt=dt)

        v.setValue(v - (dt/rho)*p.faceGrad)

    X = mesh.cellCenters[0]
    Y = mesh.cellCenters[1]
    Z = mesh.cellCenters[2]

    Vx = v[0].value
    Vy = v[1].value
    Vz = v[2].value

    fp.Viewer(vars=(v,)).plot()


    return print("Simulação completa")

if __name__ == "__main__":
    main()