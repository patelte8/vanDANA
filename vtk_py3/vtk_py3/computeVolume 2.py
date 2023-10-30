########################################################################

import sys
import numpy
import vtk

#from mat_vec_tools    import *
from vtk_py3.createFloatArray import *
from vtk_py3 import *

########################################################################


def computeVolume(pdata, orientation):

        bd = pdata.GetBounds()
        
        # Define the cutting plane
        plane=vtk.vtkPlane()
        if(orientation == 'z'):
        	plane.SetOrigin(0,0,bd[5]-0.1)
        	plane.SetNormal(0,0,-1)
        elif(orientation == 'x'):
                plane.SetOrigin(bd[0]+0.1,0,0)
                plane.SetNormal(1,0,0)
        elif(orientation == 'n'):
        	plane.SetOrigin(0,0,bd[5]+100)
        	plane.SetNormal(0,0,-1)
        
        
        # Need a plane collection for clipping
        planeCollection = vtk.vtkPlaneCollection()
        planeCollection.AddItem(plane)
        
        # The clipper generates a clipped polygonial model
        clipper = vtk.vtkClipClosedSurface()
        clipper.SetClippingPlanes(planeCollection)
        if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
                clipper.SetInput(pdata)
        else:
                clipper.SetInputData(pdata)
        clipper.SetGenerateFaces(1)
        clipper.SetScalarModeToLabels()
        clipper.Update()
        
        # Get volume using mass property
        massprop = vtk.vtkMassProperties()
        if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
                massprop.SetInput(clipper.GetOutput())
        else:
                massprop.SetInputData(clipper.GetOutput())
        return massprop.GetVolume()

if (__name__ == "__main__"):
        assert (len(sys.argv) in [2,3]), "Number of arguments must be 1 or 2. Aborting."
        pdata_mesh = readSTL(sys.argv[1])
        if (len(sys.argv) == 2):
            vol = computeVolume(pdata_mesh, "z")
        elif (len(sys.argv) == 3):
            vol = computeVolume(pdata_mesh, sys.argv[2])
        
        print("Cavity volume = ", vol)




