########################################################################

import numpy
import vtk

#from mat_vec_tools    import *
from vtk_py3.createFloatArray import *

########################################################################

def addSystolicStrains(mesh, verbose=True):

        if (verbose): print('*** addSystolicStrains ***')
        
        nb_cells = mesh.GetNumberOfCells()
        
        farray_F_dia = mesh.GetCellData().GetArray('F_dia')
        farray_F_sys = mesh.GetCellData().GetArray('F_sys')
        
        farray_E_dia = createFloatArray('E_dia', 6, nb_cells)
        farray_E_sys = createFloatArray('E_sys', 6, nb_cells)
        farray_F_num = createFloatArray('F_num', 9, nb_cells)
        farray_E_num = createFloatArray('E_num', 6, nb_cells)
        
        for num_cell in range(nb_cells):
            F_dia = numpy.reshape(farray_F_dia.GetTuple(num_cell), (3,3), order='C')
            F_sys = numpy.reshape(farray_F_sys.GetTuple(num_cell), (3,3), order='C')
            #print 'F_dia =', F_dia
            #print 'F_sys =', F_sys
        
            C = numpy.dot(numpy.transpose(F_dia), F_dia)
            E = (C - numpy.eye(3))/2
            farray_E_dia.InsertTuple(num_cell, mat_sym_to_vec_col(E))
        
            C = numpy.dot(numpy.transpose(F_sys), F_sys)
            E = (C - numpy.eye(3))/2
            farray_E_sys.InsertTuple(num_cell, mat_sym_to_vec_col(E))
        
            F = numpy.dot(F_sys, numpy.linalg.inv(F_dia))
            farray_F_num.InsertTuple(num_cell, numpy.reshape(F, 9, order='C'))
            #print 'F =', F
        
            C = numpy.dot(numpy.transpose(F), F)
            E = (C - numpy.eye(3))/2
            farray_E_num.InsertTuple(num_cell, mat_sym_to_vec_col(E))
        
        mesh.GetCellData().AddArray(farray_E_dia)
        mesh.GetCellData().AddArray(farray_E_sys)
        mesh.GetCellData().AddArray(farray_F_num)
        mesh.GetCellData().AddArray(farray_E_num)
