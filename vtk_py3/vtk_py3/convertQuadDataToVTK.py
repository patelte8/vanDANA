import vtk as vtk
from dolfin import *
import vtk_py3 as vtk_py
import math as math
import numpy as np
import os as os

def convertQuadDataToVTK(mesh, Fspace, data, filename=[]):

	nsubspace = Fspace.num_sub_spaces()
	assert (nsubspace == 3), 'Only vectorspace works'

	dim = mesh.geometry().dim()
	coord = Fspace.tabulate_dof_coordinates().reshape((-1, dim))
	npts = int(len(coord)/float(nsubspace))

	my_first, my_last = Fspace.sub(0).dofmap().ownership_range()
	
	#dofmap = [Fspace.sub(i).dofmap().dofs() for i in range(0, nsubspace)]
	x_dofs = np.arange(0, my_last-my_first, dim)
	y_dofs = np.arange(1, my_last-my_first, dim)
	z_dofs = np.arange(2, my_last-my_first, dim)
	coord_dofs = np.arange(0, (my_last-my_first)/float(dim))

	coord_reduce = coord[x_dofs]

	if(not isinstance(data, list)):
	
		points = vtk.vtkPoints()
		vec = vtk.vtkFloatArray()
		vec.SetNumberOfComponents(3)
		vec.SetName("vector")

		for x_dof, y_dof, z_dof, coord_dof in zip(x_dofs, y_dofs, z_dofs, coord_dofs):
			points.InsertNextPoint(coord_reduce[int(coord_dof)])
			#norm = math.sqrt(data.vector().array()[dofmap[0][p]]**2 +  data.vector().array()[dofmap[1][p]]**2 + data.vector().array()[dofmap[2][p]]**2)
			vec.InsertNextTuple3(data.vector().get_local()[x_dof], data.vector().get_local()[y_dof], data.vector().get_local()[z_dof])


		pdata = vtk.vtkPolyData()  
		pdata.SetPoints(points)
		pdata.GetPointData().AddArray(vec)

		glyphfilter = vtk.vtkVertexGlyphFilter()
		glyphfilter.AddInputData(pdata)
		glyphfilter.Update()

	else:

		points = vtk.vtkPoints()

		vec_array = []
		for p in range(0, len(data)):
			vec_array.append(vtk.vtkFloatArray())
			vec_array[p].SetNumberOfComponents(3)
			vec_array[p].SetName("vector"+str(p))

		for x_dof, y_dof, z_dof, coord_dof in zip(x_dofs, y_dofs, z_dofs, coord_dofs):
			points.InsertNextPoint(coord_reduce[int(coord_dof)])
			#norm = math.sqrt(data.vector().array()[dofmap[0][p]]**2 +  data.vector().array()[dofmap[1][p]]**2 + data.vector().array()[dofmap[2][p]]**2)
			for p in range(0, len(data)):
				vec_array[p].InsertNextTuple3(data[p].vector().get_local()[x_dof], \
							      data[p].vector().get_local()[y_dof], \
							      data[p].vector().get_local()[z_dof])


		pdata = vtk.vtkPolyData()  
		pdata.SetPoints(points)
		for p in range(0, len(data)):
			pdata.GetPointData().AddArray(vec_array[p])

		glyphfilter = vtk.vtkVertexGlyphFilter()
		glyphfilter.AddInputData(pdata)
		glyphfilter.Update()



	if(not (not filename)):
		filename_ = filename + str(MPI.comm_world.Get_rank()) + '.vtp'
		vtk_py.writeXMLPData(glyphfilter.GetOutput(), filename_, verbose=False)

		if(MPI.comm_world.Get_rank() == 0):
                    pvtufilename = filename + '.pvtp'
                    pvtufile = open(pvtufilename, 'w')
                    print("<?xml version=\"1.0\"?>", file=pvtufile)
                    print("<VTKFile type=\"PPolyData\" version=\"0.1\">", file=pvtufile)
                    print("<PPolyData GhostLevel=\"0\">", file=pvtufile)
                    print("<PPointData Vectors=\"vector\">", file=pvtufile)
                    print("<PDataArray type=\"Float32\" Name=\"vector\" NumberOfComponents=\"3\" />", file=pvtufile)
                    print("</PPointData>", file=pvtufile)
                    print("<PPoints>", file=pvtufile)
                    print("<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>", file=pvtufile)
                    print("</PPoints>", file=pvtufile)
                    
                    for p in range(0, MPI.comm_world.Get_size()):
                    	print("<Piece Source=\"" + os.getcwd() + "/" + filename + str(p) + '.vtp' +  "\" />", file=pvtufile)
                    print("</PPolyData>", file=pvtufile)
                    print("</VTKFile>", file=pvtufile)
                    pvtufile.close()

		
		

		
	
	
