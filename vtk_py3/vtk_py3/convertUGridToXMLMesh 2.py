########################################################################

import vtk
import dolfin
import numpy as np
from vtk.util import numpy_support

########################################################################

def convertUGridToXMLMesh(ugrid):


        num_pts = ugrid.GetNumberOfPoints()
        num_cells =  ugrid.GetNumberOfCells()
        
        celltypes = numpy_support.vtk_to_numpy(ugrid.GetCellTypesArray())
        num_tetra = np.count_nonzero(celltypes == 10)
        
        print("Number of points  = ", num_pts)
        print("Number of tetra  = ", num_tetra)
        
        mesh = dolfin.Mesh()
        editor = dolfin.MeshEditor()
        editor.open(mesh, "tetrahedron", 3, 3)  # top. and geom. dimension are both 2
        editor.init_vertices(num_pts)  # number of vertices
        editor.init_cells(num_tetra)     # number of cells
        
        for p in range(0, num_pts):
        	pt = ugrid.GetPoints().GetPoint(p)
        	editor.add_vertex(p, dolfin.Point(pt[0], pt[1], pt[2]))
        
        
        
        cnt =  0
        for p in range(0, num_cells):
        	pts = vtk.vtkIdList()
        	ugrid.GetCellPoints(p, pts)
        	if(pts.GetNumberOfIds() == 4):
        		editor.add_cell(cnt, [pts.GetId(0),  pts.GetId(1), pts.GetId(2), pts.GetId(3)])
        		cnt = cnt + 1
        	
        editor.close()
        
        return mesh
	



