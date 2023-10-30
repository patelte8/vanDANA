########################################################################

import sys
import vtk

#from mat_vec_tools import *

########################################################################

def sliceheart(domain, C, N, verbose=True):

    if (verbose): print('*** Slice Heart ***')

    plane = vtk.vtkPlane()
    plane.SetOrigin(C)
    plane.SetNormal(N)

    cutter = vtk.vtkCutter()
    cutter.SetCutFunction(plane)
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
    	cutter.SetInput(domain)
    else:
    	cutter.SetInputData(domain)
    cutter.Update()

    return cutter.GetOutput();


