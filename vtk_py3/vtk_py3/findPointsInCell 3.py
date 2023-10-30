########################################################################

import vtk

########################################################################

def findPointsInCell(points,
                         cell,
                         verbose=True):
        
        ugrid_cell = vtk.vtkUnstructuredGrid()
        ugrid_cell.SetPoints(cell.GetPoints())
        cell = vtk.vtkHexahedron()
        for k_point in range(8): cell.GetPointIds().SetId(k_point, k_point)
        cell_array_cell = vtk.vtkCellArray()
        cell_array_cell.InsertNextCell(cell)
        ugrid_cell.SetCells(vtk.VTK_HEXAHEDRON, cell_array_cell)
        
        geometry_filter = vtk.vtkGeometryFilter()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            geometry_filter.SetInputData(ugrid_cell)
        else:
            geometry_filter.SetInput(ugrid_cell)
        geometry_filter.Update()
        cell_boundary = geometry_filter.GetOutput()
        
        pdata_points = vtk.vtkPolyData()
        pdata_points.SetPoints(points)
        
        enclosed_points_filter = vtk.vtkSelectEnclosedPoints()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            enclosed_points_filter.SetSurfaceData(cell_boundary)
            enclosed_points_filter.SetInputData(pdata_points)
        else:
            enclosed_points_filter.SetSurface(cell_boundary)
            enclosed_points_filter.SetInput(pdata_points)
        enclosed_points_filter.Update()
        
        points_in_cell = [num_point for num_point in range(points.GetNumberOfPoints()) if enclosed_points_filter.GetOutput().GetPointData().GetArray('SelectedPoints').GetTuple(num_point)[0]]
        return points_in_cell
