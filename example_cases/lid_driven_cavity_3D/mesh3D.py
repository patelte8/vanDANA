from dolfin import *
from mshr import *
import numpy as np

mesh1 = BoxMesh(Point(0.0, 0.0, 0.0), Point(1.0, 1.0, 1.0), 50, 50, 50) 

# x = mesh1.coordinates()
# x[:, :2] = (x[:, :2] - 0.5) * 2
# x[:, :2] = 0.5 * (np.cos(pi * (x[:, :2] - 1.) / 2.) + 1.)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) 

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Front(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 1.0)

class Back(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)        

left = Left()
right = Right()
bottom = Bottom()
back = Back()
front = Front()
top = Top()

boundaries = MeshFunction("size_t", mesh1, mesh1.topology().dim()-1)
left.mark(boundaries, 1)
right.mark(boundaries, 3)
bottom.mark(boundaries, 2)
front.mark(boundaries, 6)
back.mark(boundaries, 5)
top.mark(boundaries, 4)

subdomains = MeshFunction("size_t", mesh1, mesh1.topology().dim())
subdomains.set_all(1)

hdf = HDF5File(mesh1.mpi_comm(), "file_f.h5", "w")
hdf.write(mesh1, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")
