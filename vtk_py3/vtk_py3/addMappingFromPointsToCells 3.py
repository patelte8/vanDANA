########################################################################

import sys
import math
import numpy
import vtk

from vtk_py3.createIntArray import *

########################################################################

def addMappingFromPointsToCells(ugrid_points, ugrid_cells, verbose=True):

        if (verbose): print('*** addMappingFromPointsToCells ***')
        
        nb_points = ugrid_points.GetNumberOfPoints()
        nb_cells = ugrid_cells.GetNumberOfCells()
        print("nb_points = " + str(nb_points))
        print("nb_cells = " + str(nb_cells))
        
        cell_locator = vtk.vtkCellLocator()
        cell_locator.SetDataSet(ugrid_cells)
        cell_locator.Update()
        
        closest_point = [0.]*3
        generic_cell = vtk.vtkGenericCell()
        num_cell = vtk.mutable(0)
        subId = vtk.mutable(0)
        dist = vtk.mutable(0.)
        
        iarray_num_cell = createIntArray("num_cell", 1, nb_points)
        
        for num_point in range(nb_points):
            point = ugrid_points.GetPoint(num_point)
        
            cell_locator.FindClosestPoint(point, closest_point, generic_cell, num_cell, subId, dist)
            #num_cell = cell_locator.FindCell(point)
        
            iarray_num_cell.InsertTuple(num_point, [num_cell])
            #print "num_point = " + str(num_point)
            #print "num_cell = " + str(num_cell)
        
        ugrid_points.GetPointData().AddArray(iarray_num_cell)
