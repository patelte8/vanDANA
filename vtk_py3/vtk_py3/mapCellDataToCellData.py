########################################################################

import numpy
import vtk

from vtk_py3.createFloatArray import *
from vtk_py3.getCellCenters   import *

########################################################################

def mapCellDataToCellData(ugrid_data,
                          ugrid_mesh,
                          field_name,
                          cell_length_ratio_for_radius=0.5,
                          threshold_min=None,
                          threshold_max=None,
                          verbose=True):

        if (verbose): print('*** mapCellDataToCellData ***')
        
        pdata_data_cell_centers = getCellCenters(ugrid_data)
        point_locator = vtk.vtkPointLocator()
        point_locator.SetDataSet(pdata_data_cell_centers)
        point_locator.Update()
        
        pdata_mesh_cell_centers = getCellCenters(ugrid_mesh)
        
        nb_components = ugrid_data.GetCellData().GetArray(field_name).GetNumberOfComponents()
        nb_cells = ugrid_mesh.GetNumberOfCells()
        farray_field_avg = createFloatArray(field_name+"_avg",
                                            nb_components,
                                            nb_cells)
        farray_field_std = createFloatArray(field_name+"_std",
                                            nb_components,
                                            nb_cells)
        
        points_within_radius = vtk.vtkIdList()
        
        for num_cell in range(nb_cells):
            l = (ugrid_mesh.GetCell(num_cell).GetLength2())**(0.5)
        
            point_locator.FindPointsWithinRadius(l*cell_length_ratio_for_radius,
                                                 pdata_mesh_cell_centers.GetPoint(num_cell),
                                                 points_within_radius)
        
            if (points_within_radius.GetNumberOfIds()):
                values = [numpy.array(ugrid_data.GetCellData().GetArray(field_name).GetTuple(points_within_radius.GetId(num_id))) for num_id in range(points_within_radius.GetNumberOfIds())]
                if (threshold_min != None):
                    values = [value for value in values if (numpy.linalg.norm(value) > threshold_min)]
                if (threshold_max != None):
                    values = [value for value in values if (numpy.linalg.norm(value) < threshold_max)]
                avg = numpy.mean(values, 0)
                std = numpy.std(values, 0)
            else:
                avg = [0]*nb_components
                std = [0]*nb_components
            farray_field_avg.InsertTuple(num_cell, avg)
            farray_field_std.InsertTuple(num_cell, std)
        
        ugrid_mesh.GetCellData().AddArray(farray_field_avg)
        ugrid_mesh.GetCellData().AddArray(farray_field_std)
