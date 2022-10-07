########################################################################

import sys
import vtk

########################################################################

def CreateVertexFromPoint(ugrid):

        vertices = vtk.vtkCellArray()
        for p in range(ugrid.GetNumberOfPoints()):
        	vert = vtk.vtkVertex()
        	vert.GetPointIds().SetId(0, p)
        	vertices.InsertNextCell(vert)
        
        ugrid.SetCells(1, vertices)



