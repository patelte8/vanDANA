from dolfin import *
from mshr import *

# Mesh
Rectangle = Rectangle(Point(0.0, 0.0), Point(2.2, 0.41))
circle = Circle(Point(0.2, 0.2), 0.05)
domain = Rectangle - circle
mesh = generate_mesh(domain, 150)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) 

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 2.2)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.41)

class Cylinder(SubDomain):
    def inside(self,x,on_boundary):
        return x[0]>0.14 and x[0]<0.26 and x[1]>0.14 and x[1]<0.26 and on_boundary

class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (-0.5, 0.5)) and between(x[0], (-1.0, 2.5)))

left = Left()
top = Top()
right = Right()
bottom = Bottom()
cylinder = Cylinder()

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundaries.set_all(0)
left.mark(boundaries, 1)
top.mark(boundaries, 4)
right.mark(boundaries, 3)
bottom.mark(boundaries, 2)
cylinder.mark(boundaries, 5)

obstacle = Obstacle()

subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
subdomains.set_all(0)
obstacle.mark(subdomains, 1)

x = mesh.coordinates()
Lsc = 0.1
x[:, :] /= Lsc

hdf = HDF5File(mesh.mpi_comm(), "file_f.h5", "w")
hdf.write(mesh, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")
