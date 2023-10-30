########################################################################

import numpy
import vtk

from vtk_py3.createFloatArray import *

########################################################################

def mapPointDataToCellData(ugrid_data,
                           ugrid_mesh,
                           field_name,
                           cell_length_ratio_for_radius=0.5,
                           verbose=True):

        if (verbose): print('*** mapPointDataToCellData ***')
        
        nb_cells = ugrid_mesh.GetNumberOfCells()
        
        filter_cell_centers = vtk.vtkCellCenters()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            filter_cell_centers.SetInputData(ugrid_mesh)
        else:
            filter_cell_centers.SetInput(ugrid_mesh)
        filter_cell_centers.Update()
        pdata_mesh_cell_centers = filter_cell_centers.GetOutput()
        
        point_locator = vtk.vtkPointLocator()
        point_locator.SetDataSet(ugrid_data)
        point_locator.Update()
        
        nb_components = ugrid_data.GetPointData().GetArray(field_name).GetNumberOfComponents()
        farray_tensors = createFloatArray(field_name,
                                          nb_components,
                                          nb_cells)
        
        points_within_radius = vtk.vtkIdList()
        
        for num_cell in range(nb_cells):
            l = (ugrid_mesh.GetCell(num_cell).GetLength2())**(0.5)
        
            point_locator.FindPointsWithinRadius(l*cell_length_ratio_for_radius,
                                                 pdata_mesh_cell_centers.GetPoint(num_cell),
                                                 points_within_radius)
            #points_in_cell = findPointsInCell(ugrid_data.GetPoints(), ugrid_mesh.GetCell(num_cell))
        
            if (points_within_radius.GetNumberOfIds()):
                tensor = sum([numpy.array(ugrid_data.GetPointData().GetArray(field_name).GetTuple(points_within_radius.GetId(num_id))) for num_id in range(points_within_radius.GetNumberOfIds())])/points_within_radius.GetNumberOfIds()
            else:
                tensor = [0]*nb_components
            farray_tensors.InsertTuple(num_cell, tensor)
        
        ugrid_mesh.GetCellData().AddArray(farray_tensors)
