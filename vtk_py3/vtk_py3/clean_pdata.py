########################################################################

import sys
import vtk

#from mat_vec_tools import *

########################################################################

def clean_pdata(pdata):

        cleanpdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        	cleanpdata.SetInputData(pdata)
        else:
        	cleanpdata.SetInput(pdata)
        cleanpdata.Update()
        
        return cleanpdata.GetOutput()
	


