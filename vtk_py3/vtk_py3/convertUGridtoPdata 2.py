########################################################################

import vtk

########################################################################


def convertUGridtoPdata(ugrid):
    
        geometry = vtk.vtkGeometryFilter()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
                geometry.SetInputData(ugrid)
        else:
                geometry.SetInput(ugrid)

        geometry.Update()

        return geometry.GetOutput()

