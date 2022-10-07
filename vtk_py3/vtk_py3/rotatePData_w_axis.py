########################################################################

import numpy
import vtk
from pyquaternion import Quaternion

########################################################################

def rotatePData_w_axis(old_pdata, angle, axis, verbose=True):

    if (verbose): print('*** rotatePData_w_axis ***')

    nb_points = old_pdata.GetNumberOfPoints()

    new_points = vtk.vtkPoints()
    new_points.SetNumberOfPoints(nb_points)

    my_quaternion = Quaternion(axis=axis, angle=angle)

    old_point = numpy.array([0.]*3)

    for num_point in range(nb_points):
        old_pdata.GetPoint(num_point, old_point)
        #print old_point

        new_point = my_quaternion.rotate(old_point)
        #print new_point

        new_points.InsertPoint(num_point, new_point)

    new_pdata = vtk.vtkPolyData()
    new_pdata.SetPoints(new_points)
    new_pdata.SetPolys(old_pdata.GetPolys())
    new_pdata.SetLines(old_pdata.GetLines())
    
    return new_pdata
 
