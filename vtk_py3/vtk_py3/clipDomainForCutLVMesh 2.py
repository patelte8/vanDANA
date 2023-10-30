########################################################################

import sys
import vtk

#from mat_vec_tools import *

########################################################################

def clipDomainForCutLVMesh(domain, C, N, verbose=True):

        if (verbose): print('*** clipDomainForCutLVMesh ***')
        
        plane = vtk.vtkPlane()
        plane.SetOrigin(C)
        plane.SetNormal(N)
        
        clip = vtk.vtkClipPolyData()
        clip.SetClipFunction(plane)
        clip.GenerateClippedOutputOn()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            clip.SetInputData(domain)
        else:
            clip.SetInput(domain)
        clip.Update()
        clipped0 = clip.GetOutput(0)
        clipped1 = clip.GetOutput(1)
        
        if (clipped0.GetNumberOfPoints() > clipped1.GetNumberOfPoints()):
            return clipped0
        else:
            return clipped1
