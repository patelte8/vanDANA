########################################################################

import vtk

########################################################################

def getPDataNormals(pdata,
                    flip=False,
                    verbose=True):

    if (verbose): print('*** getPDataNormals ***')

    poly_data_normals = vtk.vtkPolyDataNormals()
    poly_data_normals.ComputePointNormalsOff()
    poly_data_normals.ComputeCellNormalsOn()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        poly_data_normals.SetInputData(pdata)
    else:
        poly_data_normals.SetInput(pdata)
    if (flip): poly_data_normals.FlipNormalsOn()
    else:      poly_data_normals.FlipNormalsOff()
    poly_data_normals.Update()

    return poly_data_normals.GetOutput()
