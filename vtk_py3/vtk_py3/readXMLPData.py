########################################################################

import vtk

########################################################################

def readXMLPData(pdata_file_name, verbose=True):
        if (verbose): print('*** readXMLPData ***')
        
        pdata_reader = vtk.vtkXMLPolyDataReader()
        pdata_reader.SetFileName(pdata_file_name)
        pdata_reader.Update()
        pdata = pdata_reader.GetOutput()
        
        if (verbose):
            print("nb_points = " + str(pdata.GetNumberOfPoints()))
            print("nb_verts = " + str(pdata.GetNumberOfVerts()))
            print("nb_lines = " + str(pdata.GetNumberOfLines()))
            print("nb_polys = " + str(pdata.GetNumberOfPolys()))
            print("nb_strips = " + str(pdata.GetNumberOfStrips()))
        
        return pdata
