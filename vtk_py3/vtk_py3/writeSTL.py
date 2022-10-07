########################################################################

import vtk

########################################################################

def writeSTL(stl, stl_file_name, verbose=True):
    if (verbose): print('*** writeSTL ***')

    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetFileName(stl_file_name)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        stl_writer.SetInputData(stl)
    else:
        stl_writer.SetInput(stl)
    stl_writer.Update()
    stl_writer.Write()
