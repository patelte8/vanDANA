import dolfin as dolfin
import vtk_py3 as vtk_py
import vtk

casename = "ellipsoidal"
meshname = "ellipsoidal"
directory = "./"
meshsize = 0.05
isepiflip = False
isendoflip = True

# Test create_ellipsoidal_LV
ugrid = vtk_py.create_ellipsoidal_LV(gmshcmd="gmsh", meshsize=meshsize)


xmlgrid = vtk_py.convertUGridToXMLMesh(ugrid)

xmlgrid, xmlfacet, xmledges = vtk_py.extractFeNiCsBiVFacet(ugrid, geometry="LV")

VQuadelem = dolfin.VectorElement("Quadrature", 
                              xmlgrid.ufl_cell(), 
                              degree=4, 
                              quad_scheme="default")
VQuadelem._quad_scheme = 'default'

fiberFS = dolfin.FunctionSpace(xmlgrid, VQuadelem)

ef, es, en, eC, eL, eR = vtk_py.addLVfiber(xmlgrid, fiberFS, casename, 60, -60, [] , isepiflip, isendoflip)

f = dolfin.HDF5File(xmlgrid.mpi_comm(), directory + meshname+".hdf5", 'w')
f.write(xmlgrid, casename)
f.close()

f = dolfin.HDF5File(xmlgrid.mpi_comm(), directory + meshname+".hdf5", 'a') 
f.write(xmlfacet, casename+"/"+"facetboundaries") 
f.write(xmledges, casename+"/"+"edgeboundaries") 
f.write(ef, casename+"/"+"eF") 
f.write(es, casename+"/"+"eS") 
f.write(en, casename+"/"+"eN")
f.write(eC, casename+"/"+"eC") 
f.write(eL, casename+"/"+"eL") 
f.write(eR, casename+"/"+"eR") 
f.close()

dolfin.File("facetboundaries.pvd")  << xmlfacet
dolfin.File("Edges.pvd")  << xmledges

