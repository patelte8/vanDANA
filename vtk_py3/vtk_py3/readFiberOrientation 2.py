########################################################################

import numpy
import vtk

from vtk_py3.createFloatArray import *

########################################################################

def readFiberOrientation(mesh,
                     filename,
                     verbose=True):

        if (verbose): print('*** readFiberOrientation ***')
        
        nb_cells = mesh.GetNumberOfCells()
        
        if verbose: print('Reading fiber orientations...')
        
        eF_array = createFloatArray('eF', 3, nb_cells)
        eS_array = createFloatArray('eS', 3, nb_cells)
        eN_array = createFloatArray('eN', 3, nb_cells)
        
        file = open(filename, 'r')
        file.readline()
        
        num_cell = 0
        for line in file:
            line = line.split(', ')
            #print line
        
            eF = [float(item) for item in line[1:4]]
            eS = [float(item) for item in line[4:7]]
            eN = numpy.cross(eF,eS)
            #print "eF =", eF
            #print "eS =", eS
            #print "eN =", eN
        
            eF_array.InsertTuple(num_cell, eF)
            eS_array.InsertTuple(num_cell, eS)
            eN_array.InsertTuple(num_cell, eN)
        
            num_cell += 1
            #print "num_cell =", num_cell
        
        file.close()
        
        mesh.GetCellData().AddArray(eF_array)
        mesh.GetCellData().AddArray(eS_array)
        mesh.GetCellData().AddArray(eN_array)
        
        return mesh
