import vtk as vtk
from dolfin import *
import vtk_py as vtk_py
import math as math
import numpy as np
import os as os

def convertQuadScalarDataToVTK(mesh, Fspace, data, filename=[]):

	nsubspace = Fspace.num_sub_spaces()
	assert (nsubspace == 0), 'Only scalar space works'

	dim = mesh.geometry().dim()
	coord = Fspace.tabulate_dof_coordinates().reshape((-1, dim))
	npts = int(len(coord))

	my_first, my_last = Fspace.dofmap().ownership_range()
	
	x_dofs = np.arange(0, my_last-my_first)
	coord_dofs = np.arange(0, (my_last-my_first))

	coord_reduce = coord[x_dofs]

	if(not isinstance(data, list)):
	
		points = vtk.vtkPoints()
		scalar = vtk.vtkFloatArray()
		scalar.SetNumberOfComponents(1)
		scalar.SetName("scalar")

		for x_dof, coord_dof in zip(x_dofs, coord_dofs):
			points.InsertNextPoint(coord_reduce[int(coord_dof)])
			scalar.InsertNextValue(data.vector().array()[x_dof])


		pdata = vtk.vtkPolyData()  
		pdata.SetPoints(points)
		pdata.GetPointData().AddArray(scalar)

		glyphfilter = vtk.vtkVertexGlyphFilter()
		glyphfilter.AddInputData(pdata)
		glyphfilter.Update()

	else:

		points = vtk.vtkPoints()

		scalar_array = []
		for p in range(0, len(data)):
			scalar_array.append(vtk.vtkFloatArray())
			scalar_array[p].SetNumberOfComponents(3)
			scalar_array[p].SetName("scalar"+str(p))

		for x_dof, coord_dof in zip(x_dofs, coord_dofs):
			points.InsertNextPoint(coord_reduce[int(coord_dof)])
			for p in range(0, len(data)):
				scalar_array[p].InsertNextValue(data[p].vector().array()[x_dof])


		pdata = vtk.vtkPolyData()  
		pdata.SetPoints(points)
		for p in range(0, len(data)):
			pdata.GetPointData().AddArray(scalar_array[p])

		glyphfilter = vtk.vtkVertexGlyphFilter()
		glyphfilter.AddInputData(pdata)
		glyphfilter.Update()



	if(not (not filename)):
		filename_ = filename + str(MPI.rank(mpi_comm_world())) + '.vtp'
		#vtk_py.writeXMLPData(glyphfilter.GetOutput(), filename_, verbose=False)
		vtk_py.writeXMLPData(pdata, filename_, verbose=False)

		if(MPI.rank(mpi_comm_world()) == 0):
			pvtufilename = filename + '.pvtp'
			pvtufile = open(pvtufilename, 'w')
			print >>pvtufile, "<?xml version=\"1.0\"?>"
			print >>pvtufile, "<VTKFile type=\"PPolyData\" version=\"0.1\">"
		  	print >>pvtufile, "<PPolyData GhostLevel=\"0\">"
			print >>pvtufile, "<PPointData Vectors=\"vector\">"
		  	print >>pvtufile, "<PDataArray type=\"Float32\" Name=\"vector\" NumberOfComponents=\"3\" />"
			print >>pvtufile, "</PPointData>"
			print >>pvtufile, "<PPoints>"
			print >>pvtufile, "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"
			print >>pvtufile, "</PPoints>"
		
			for p in range(0, MPI.size(mpi_comm_world())):
		    		print >>pvtufile, "<Piece Source=\"" + os.getcwd() + "/" + filename + str(p) + '.vtp' +  "\" />"
		  	print >>pvtufile, "</PPolyData>"
			print >>pvtufile, "</VTKFile>"
			pvtufile.close()

		
		

	return pdata		
	
	

