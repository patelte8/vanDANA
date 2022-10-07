from dolfin import *
from mshr import *

mesh1 = RectangleMesh(Point(0.0, 0.0), Point(11, 4.1), 700, 260) 


class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) 

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 11.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 4.1)


class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (-1.0, 5)) and between(x[0], (-0.1, 12)))

left = Left()
top = Top()
right = Right()
bottom = Bottom()

boundaries = MeshFunction("size_t", mesh1, mesh1.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 4)
right.mark(boundaries, 3)
bottom.mark(boundaries, 2)

obstacle = Obstacle()

subdomains = MeshFunction("size_t", mesh1, mesh1.topology().dim())
subdomains.set_all(0)
obstacle.mark(subdomains, 1)

hdf = HDF5File(mesh1.mpi_comm(), "file_f.h5", "w")
hdf.write(mesh1, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")
