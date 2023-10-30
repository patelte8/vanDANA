#coding=utf8

########################################################################

import numpy
import vtk

########################################################################

def getABPointsFromTTTSectors(ugrid_sectors,
                                  verbose=True):
        
        if (verbose): print('*** getABPointsFromTTTSectors ***')
        
        nb_points = ugrid_sectors.GetNumberOfPoints()
        nb_cells = ugrid_sectors.GetNumberOfCells()
        
        nb_csects = 12
        nb_rsects = 3
        nb_slices = nb_points / (nb_rsects+1) / (nb_csects+1)
        if (verbose): print("nb_slices =", nb_slices)
        
        zmin = ugrid_sectors.GetPoint(0)[2]
        zmax = ugrid_sectors.GetPoint(ugrid_sectors.GetNumberOfPoints()-1)[2]
        
        dist_btw_slices = abs(zmin-zmax)/(nb_slices-1)
        if (verbose): print("dist_btw_slices =", dist_btw_slices)
        
        A = ugrid_sectors.GetPoints().GetPoint(0)
        B = ugrid_sectors.GetPoints().GetPoint(6)
        C = ugrid_sectors.GetPoints().GetPoint(3)
        D = ugrid_sectors.GetPoints().GetPoint(9)
        
        #print A
        #print B
        #print C
        #print D
        
        Px = ((A[0]*B[1]-A[1]*B[0])*(C[0]-D[0])-(A[0]-B[0])*(C[0]*D[1]-C[1]*D[0]))/((A[0]-B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]-D[0]))
        Py = ((A[0]*B[1]-A[1]*B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]*D[1]-C[1]*D[0]))/((A[0]-B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]-D[0]))
        
        #print Px
        #print Py
        
        A = [Px, Py, zmin]
        B = [Px, Py, zmax]
        
        if (verbose): print("A =", A)
        if (verbose): print("B =", B)
        
        points_AB = vtk.vtkPoints()
        points_AB.InsertNextPoint(A)
        points_AB.InsertNextPoint(B)
        
        cells_AB = vtk.vtkCellArray()
        cell_AB  = vtk.vtkVertex()
        cell_AB.GetPointIds().SetId(0, 0)
        cells_AB.InsertNextCell(cell_AB)
        cell_AB.GetPointIds().SetId(0, 1)
        cells_AB.InsertNextCell(cell_AB)
        
        AB_ugrid = vtk.vtkUnstructuredGrid()
        AB_ugrid.SetPoints(points_AB)
        AB_ugrid.SetCells(vtk.VTK_VERTEX, cells_AB)
        
        return points_AB



