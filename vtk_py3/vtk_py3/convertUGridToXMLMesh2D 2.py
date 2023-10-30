########################################################################

import vtk
import dolfin
import numpy as np

########################################################################

def convertUGridToXMLMesh2D(ugrid):


        num_pts = ugrid.GetNumberOfPoints()
	num_cells =  ugrid.GetNumberOfCells()

	print num_pts
	print num_cells

	celltype = vtk.vtkCellTypes()
	ugrid.GetCellTypes(celltype)
	celltype =  celltype.GetCellType(0)

        mesh = dolfin.Mesh()
	editor = dolfin.MeshEditor()
	editor.open(mesh, 2, 2)  # top. and geom. dimension are both 2
	editor.init_vertices(num_pts)  # number of vertices

	if(celltype == 8):
		editor.init_cells(2*num_cells)     # number of cells
	else:	
		editor.init_cells(num_cells)     # number of cells

	for p in range(0, num_pts):
		pt = ugrid.GetPoints().GetPoint(p)
		editor.add_vertex(p, pt[0], pt[1])


	for p in range(0, num_cells):
		pts = vtk.vtkIdList()
		ugrid.GetCellPoints(p, pts)

		if(celltype == 8):
			editor.add_cell(2*p, pts.GetId(0),  pts.GetId(1), pts.GetId(2))
			editor.add_cell(2*p+1, pts.GetId(1),  pts.GetId(2), pts.GetId(3))
		else:
			editor.add_cell(p, pts.GetId(0),  pts.GetId(1), pts.GetId(2))
		
	editor.close()

	return mesh
	



