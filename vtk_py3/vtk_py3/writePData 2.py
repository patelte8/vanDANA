########################################################################

import vtk

########################################################################

def writePData(pdata, pdata_file_name, verbose=True):
    if (verbose): print('*** writePData ***')

    pdata_writer = vtk.vtkPolyDataWriter()
    pdata_writer.SetFileName(pdata_file_name)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        pdata_writer.SetInputData(pdata)
    else:
        pdata_writer.SetInput(pdata)
    pdata_writer.Update()
    pdata_writer.Write()
