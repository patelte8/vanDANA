########################################################################

import vtk

########################################################################

def createFloatArray(name, nb_components=1, nb_tuples=0, verbose=True):
    farray = vtk.vtkFloatArray()
    farray.SetName(name)
    farray.SetNumberOfComponents(nb_components)
    farray.SetNumberOfTuples(nb_tuples)
    return farray