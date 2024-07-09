from dolfin import *
from mshr import *

# Mesh
Rectangle = Rectangle(Point(-10.0, -15.0), Point(35.0, 15.0))
circle = Circle(Point(0.0, 0.0), 0.5)
domain = Rectangle - circle
mesh1 = generate_mesh(domain, 500)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], -10.0) 

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 35.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], -15.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 15.0)

class Cylinder(SubDomain):
    def inside(self,x,on_boundary):
        return x[0]>-0.6 and x[0]<0.6 and x[1]>-0.6 and x[1]<0.6 and on_boundary

class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (-20.0, 20)) and between(x[0], (-15.0, 40)))

left = Left()
top = Top()
right = Right()
bottom = Bottom()
cylinder = Cylinder()

boundaries = MeshFunction("size_t", mesh1, mesh1.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 4)
right.mark(boundaries, 3)
bottom.mark(boundaries, 2)
cylinder.mark(boundaries, 5)

obstacle = Obstacle()

subdomains = MeshFunction("size_t", mesh1, mesh1.topology().dim())
subdomains.set_all(0)
obstacle.mark(subdomains, 1)

hdf = HDF5File(mesh1.mpi_comm(), "file_f.h5", "w")
hdf.write(mesh1, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")