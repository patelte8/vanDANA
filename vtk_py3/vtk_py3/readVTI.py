########################################################################

import vtk

########################################################################

def readVTI(vti_file_name, verbose=True):
    if (verbose): print('Reading VTK XML mesh...')

    image = vtk.vtkImageData()
    image_reader = vtk.vtkXMLImageDataReader()
    image_reader.SetFileName(vti_file_name)
    image_reader.Update()
    
    image.DeepCopy(image_reader.GetOutput())

    if (verbose):
        nb_cells = image.GetNumberOfCells()
        print('nb_cells =', nb_cells)

    return image

