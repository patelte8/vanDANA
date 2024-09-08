from dolfin import *
from mshr import *

mesh1 = BoxMesh(Point(0.0, -1., -1.), Point(5.0, 1., 1.), 60, 24, 24)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return (x[1]*x[1] + x[2]*x[2] >= 0.25) and near(x[0], 0.0) 

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return (x[1]*x[1] + x[2]*x[2] >= 0.25) and near(x[0], 5.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], -1.)

class Front(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 1.)

class Back(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], -1.)

class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return (x[1]*x[1] + x[2]*x[2] < 0.25) and near(x[0], 5.0)

class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return (x[1]*x[1] + x[2]*x[2] < 0.25) and near(x[0], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.)        

left = Left()
right = Right()
bottom = Bottom()
back = Back()
inlet = Inlet()
outlet = Outlet()
front = Front()
top = Top()

boundaries = MeshFunction("size_t", mesh1, mesh1.topology().dim()-1)
inlet.mark(boundaries, 1)
outlet.mark(boundaries, 2)
left.mark(boundaries, 3)
right.mark(boundaries, 4)
bottom.mark(boundaries, 5)
front.mark(boundaries, 6)
back.mark(boundaries, 7)
top.mark(boundaries, 8)

# file = File("Mesh_f.pvd") 
# file << mesh1
# file << boundaries

class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (-1.7, 1.5)) and between(x[0], (-1., 5.1)) and between(x[2], (-1.5, 1.5)))


subdomains = MeshFunction("size_t", mesh1, mesh1.topology().dim())
subdomains.set_all(0)
obstacle = Obstacle()
obstacle.mark(subdomains, 1)

hdf = HDF5File(mesh1.mpi_comm(), "file_f.h5", "w")
hdf.write(mesh1, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")