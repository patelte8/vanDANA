########################################################################

import numpy
import vtk

from createFloatArray import *

########################################################################

def readDynaDeformationGradients(mesh,
                                 hystory_files_basename,
                                 array_name,
                                 verbose=True):

    if (verbose): print '*** readDynaDeformationGradients ***'

    nb_cells = mesh.GetNumberOfCells()

    history_files_names = [hystory_files_basename + '.history#' + str(num) for num in range(11,20)]

    F_list = [[0. for num_component in range(9)] for num_cell in range(nb_cells)]

    for num_component in range(9):
        history_file = open(history_files_names[num_component], 'r')
        for line in history_file:
            if line.startswith('*') or line.startswith('$'): continue
            line = line.split()
            F_list[int(line[0])-1][num_component] = float(line[1])
        history_file.close()

    F_array = createFloatArray(array_name, 9, nb_cells)

    for num_cell in range(nb_cells):
        F_array.InsertTuple(num_cell, F_list[num_cell])

    if (verbose): print "nb_tuples = " + str(F_array.GetNumberOfTuples())

    mesh.GetCellData().AddArray(F_array)
