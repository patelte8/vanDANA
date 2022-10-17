########################################################################

import vtk

########################################################################

def readXMLPUGrid(mesh_file_name, verbose=True):
    if (verbose): print('*** readXMLPUGrid ***')

    mesh_reader = vtk.vtkXMLPUnstructuredGridReader()
    mesh_reader.SetFileName(mesh_file_name)
    mesh_reader.Update()
    mesh = mesh_reader.GetOutput()

    if (verbose):
        nb_points = mesh.GetNumberOfPoints()
        print('nb_points =', nb_points)

        nb_cells = mesh.GetNumberOfCells()
        print('nb_cells =', nb_cells)

    return mesh
