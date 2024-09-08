from dolfin import *

mesh2 = Mesh("pipe.xml")
subdomains2 = MeshFunction("size_t", mesh2, "pipe_physical_region.xml")
boundaries2 = MeshFunction("size_t", mesh2, "pipe_facet_region.xml")

hdf = HDF5File(mesh2.mpi_comm(), "file_s.h5", "w")
hdf.write(mesh2, "/mesh")
hdf.write(subdomains2, "/subdomains")
hdf.write(boundaries2, "/boundaries")