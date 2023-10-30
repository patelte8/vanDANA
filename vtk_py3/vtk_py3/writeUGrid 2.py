########################################################################

import vtk

########################################################################

def writeUGrid(ugrid, ugrid_file_name, verbose=True):
    if (verbose): print('*** writeUGrid ***')

    ugrid_writer = vtk.vtkUnstructuredGridWriter()
    ugrid_writer.SetFileName(ugrid_file_name)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        ugrid_writer.SetInputData(ugrid)
    else:
        ugrid_writer.SetInput(ugrid)
    ugrid_writer.Update()
    ugrid_writer.Write()
