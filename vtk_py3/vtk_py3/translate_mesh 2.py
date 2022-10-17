########################################################################

import sys
import numpy
import vtk
from vtk_py3 import *

########################################################################

def translate_mesh(ugrid, vec):

        transform = vtk.vtkTransform()
        transform.Translate(vec[0],vec[1],vec[2])
        transform.Update()
        
        transformfilter = vtk.vtkTransformFilter()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        	transformfilter.SetInputData(ugrid)
        else:
        	transformfilter.SetInput(ugrid)
        transformfilter.SetTransform(transform)
        transformfilter.Update()
        
        return transformfilter.GetOutput()


 



