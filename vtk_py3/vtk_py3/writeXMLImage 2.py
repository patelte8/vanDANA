########################################################################

import vtk

########################################################################

def writeXMLImage(image, image_file_name, verbose=True):
    if (verbose): print('*** writeXMLImage ***')

    image_writer = vtk.vtkXMLImageDataWriter()
    image_writer.SetFileName(image_file_name)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        image_writer.SetInputData(image)
    else:
        image_writer.SetInput(image)
    image_writer.Update()
    image_writer.Write()
