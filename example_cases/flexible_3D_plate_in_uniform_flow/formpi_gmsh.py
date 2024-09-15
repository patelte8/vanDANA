from dolfin import *
from mshr import *

mesh2 = Mesh("flag3D.xml")
subdomains2 = MeshFunction("size_t", mesh2, "flag3D_physical_region.xml")
boundaries2 = MeshFunction("size_t", mesh2, "flag3D_facet_region.xml")

hdf = HDF5File(mesh2.mpi_comm(), "file_s.h5", "w")
hdf.write(mesh2, "/mesh")
hdf.write(subdomains2, "/subdomains")
hdf.write(boundaries2, "/boundaries")