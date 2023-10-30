########################################################################

import sys
import vtk

#from mat_vec_tools import *

########################################################################

def clipheart(domain, C, N, isinsideout, verbose=True):

    if (verbose): print('*** Slice Heart ***')

    plane = vtk.vtkPlane()
    plane.SetOrigin(C)
    plane.SetNormal(N)

    clipper = vtk.vtkClipPolyData()
    clipper.SetClipFunction(plane)
    if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
    	clipper.SetInput(domain)
    else:
    	clipper.SetInputData(domain)
    clipper.SetInsideOut(isinsideout)
    clipper.Update()

    return clipper.GetOutput();


