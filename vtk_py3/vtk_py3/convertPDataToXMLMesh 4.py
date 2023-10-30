########################################################################

import vtk
import dolfin
import numpy as np
from vtk.util import numpy_support

########################################################################

def convertPDataToXMLMesh(pdata):


        num_pts = pdata.GetNumberOfPoints()
        num_cells =  pdata.GetNumberOfCells()
        
        #celltypes = numpy_support.vtk_to_numpy(pdata.GetCellTypesArray())
        #num_triangle = np.count_nonzero(celltypes == 5)
        num_triangle = num_cells
        
        print("Number of points  = ", num_pts)
        print("Number of triangle  = ", num_triangle)
        
        mesh = dolfin.Mesh()
        editor = dolfin.MeshEditor()
        editor.open(mesh, 'triangle', 2, 3)  # top. and geom. dimension are both 2
        editor.init_vertices(num_pts)  # number of vertices
        editor.init_cells(num_triangle)     # number of cells
        
        for p in range(0, num_pts):
        	pt = pdata.GetPoints().GetPoint(p)
        	editor.add_vertex(p, dolfin.Point(pt[0], pt[1], pt[2]))
        
        
        cnt =  0
        for p in range(0, num_cells):
        	pts = vtk.vtkIdList()
        	pdata.GetCellPoints(p, pts)
        	if(pts.GetNumberOfIds() == 3):
        		editor.add_cell(cnt, [pts.GetId(0),  pts.GetId(1), pts.GetId(2)])
        		cnt = cnt + 1
        	
        editor.close()
        
        return mesh
	



