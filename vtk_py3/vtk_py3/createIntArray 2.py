########################################################################

import vtk

########################################################################

def createIntArray(name, nb_components, nb_tuples, verbose=True):
    iarray = vtk.vtkIntArray()
    iarray.SetName(name)
    iarray.SetNumberOfComponents(nb_components)
    iarray.SetNumberOfTuples(nb_tuples)
    return iarray
