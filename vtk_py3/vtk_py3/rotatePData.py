########################################################################

import numpy
import vtk

########################################################################

def rotatePData(old_pdata, C, R, verbose=True):
    if (verbose): print '*** rotatePData ***'

    nb_points = old_pdata.GetNumberOfPoints()

    new_points = vtk.vtkPoints()
    new_points.SetNumberOfPoints(nb_points)

    old_point = numpy.array([0.]*3)

    for num_point in range(nb_points):
        old_pdata.GetPoint(num_point, old_point)
        #print old_point

        new_point = C + numpy.dot(R, old_point - C)
        #print new_point

        new_points.InsertPoint(num_point, new_point)

    new_pdata = vtk.vtkPolyData()
    new_pdata.SetPoints(new_points)
    new_pdata.SetPolys(old_pdata.GetPolys())
    
    return new_pdata
