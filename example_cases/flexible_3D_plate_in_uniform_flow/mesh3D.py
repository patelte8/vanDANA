from dolfin import *
from mshr import *

mesh1 = BoxMesh(Point(0.0, 0.0, 0.0), Point(4.0, 2., 1.), 60, 30, 15)

# ---------------------------------------------------------------------------
mesh2 = Mesh()
hdf2 = HDF5File(mesh2.mpi_comm(), "file_s.h5", "r")
hdf2.read(mesh2, "/mesh", False)

bmesh = BoundaryMesh(mesh2, 'exterior')

tree = bmesh.bounding_box_tree()

# ---------------------------------------------------------------------------

V = FunctionSpace(mesh1, 'CG', 1)
v_2_d = vertex_to_dof_map(V)
bdry_distance = Function(V)

nodes_initial  = MPI.sum(mesh1.mpi_comm(), bdry_distance.vector().get_local().size)
print(nodes_initial)

values = bdry_distance.vector().get_local()

for index, vertex in enumerate(vertices(mesh1)):
    _, d = tree.compute_closest_entity(vertex.point())
    values[v_2_d[index]] = d

bdry_distance.vector().set_local(values)
bdry_distance.vector().apply('insert')
bdry_distance.set_allow_extrapolation(True)

# Number of refinements
dist = 0.35
nor = 2

for i in range(1, nor+1):

    # Selecting edges to refine
    class Border(SubDomain):
        def inside(self, x, on_boundary):
            #return True if bdry_distance(x) <= (dist/i) else False # or (x[0]<5. and x[0]>-1. and x[1]<2. and x[1]>-2. and x[2]<0. and x[2]>-8.)
            return True if (x[0]<2.3 and x[0]>0.8 and x[1]<1.5 and x[1]>0.5 and x[2]<1.1 and x[2]>-0.1) else False
    Border = Border()
    edge_markers = MeshFunction("bool", mesh1, 1, False)
    Border.mark(edge_markers, True)

    mesh1 = refine(mesh1, edge_markers)

#mesh1.smooth(20)

V = FunctionSpace(mesh1, 'CG', 1)
vector_new = Function(V)
nodes_final  = MPI.sum(mesh1.mpi_comm(), vector_new.vector().get_local().size)
print(nodes_final)

# ---------------------------------------------------------------------------

class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0) 

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 4.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1],0.0)

class Front(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 1.0)

class Back(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[2], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 2.0)        

left = Left()
right = Right()
bottom = Bottom()
back = Back()
front = Front()
top = Top()

boundaries = MeshFunction("size_t", mesh1, mesh1.topology().dim()-1)
left.mark(boundaries, 1)
right.mark(boundaries, 4)
bottom.mark(boundaries, 2)
front.mark(boundaries, 3)
back.mark(boundaries, 6)
top.mark(boundaries, 5)

class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (between(x[1], (-1.7, 2.5)) and between(x[0], (-1.1, 5.1)) and between(x[2], (-1.1, 1.1)))


subdomains = MeshFunction("size_t", mesh1, mesh1.topology().dim())
subdomains.set_all(0)
obstacle = Obstacle()
obstacle.mark(subdomains, 1)

# -----------------------------------------------------------------------

# file = File("Mesh_f.pvd") 
# file << mesh1
# file << boundaries

hdf = HDF5File(mesh1.mpi_comm(), "file_f.h5", "w")
hdf.write(mesh1, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")