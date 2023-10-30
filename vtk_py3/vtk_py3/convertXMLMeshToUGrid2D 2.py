########################################################################

import vtk
import dolfin

########################################################################

def convertXMLMeshToUGrid2D(mesh):

	connectivity =  mesh.cells()
	coords =  mesh.coordinates()

	points = vtk.vtkPoints()
	for coord in coords:
		points.InsertNextPoint(coord[0], coord[1], coord[2])
	
	cellarray = vtk.vtkCellArray()
	for cell in connectivity:
		triangle = vtk.vtkTriangle()
		triangle.GetPointIds().SetId(0, cell[0])
		triangle.GetPointIds().SetId(1, cell[1])
		triangle.GetPointIds().SetId(2, cell[2])

		cellarray.InsertNextCell(triangle)

	ugrid = vtk.vtkUnstructuredGrid()
	ugrid.SetPoints(points)
	ugrid.SetCells(triangle.GetCellType(), cellarray)

	return ugrid	
	



