########################################################################

import vtk

########################################################################

def readXMLImage(image_file_name, verbose=True):
        if (verbose): print('*** readXMLImage ***')
        
        image_reader = vtk.vtkXMLImageDataReader()
        image_reader.SetFileName(image_file_name)
        image_reader.Update()
        image = image_reader.GetOutput()
        
        if (verbose):
            print("nb_points = " + str(image.GetNumberOfPoints()))
        
        return image
