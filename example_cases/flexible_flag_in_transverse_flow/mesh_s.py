from dolfin import *
from mshr import *

mesh2 = RectangleMesh(Point(3.9894, 0.0), Point(4.0106, 0.8), 3, 80) 


class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (-0.2, 1.0)) and between(x[0], (3.0, 5.0)))

bottom = Bottom()

boundaries = MeshFunction("size_t", mesh2, mesh2.topology().dim()-1)
boundaries.set_all(0)
bottom.mark(boundaries, 1)

obstacle = Obstacle()

subdomains = MeshFunction("size_t", mesh2, mesh2.topology().dim())
subdomains.set_all(0)
obstacle.mark(subdomains, 1)

hdf = HDF5File(mesh2.mpi_comm(), "file_s.h5", "w")
hdf.write(mesh2, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")
