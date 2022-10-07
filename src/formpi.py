from dolfin import *

# Create mesh
mesh = Mesh("LA_CB.xml")
boundaries = MeshFunction("size_t", mesh, "LA_CB_facet_region.xml") 

x = mesh.coordinates()
scaling_factor = 0.0001
x[:, :] *= scaling_factor

hdf = HDF5File(mesh.mpi_comm(), "file.h5", "w")
hdf.write(mesh, "/mesh")
hdf.write(boundaries, "/boundaries")

