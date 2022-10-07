########################################################################

import numpy
import vtk

#from mat_vec_tools    import *
from vtk_py3.createFloatArray import *

########################################################################

def rotateSymmetricMatrix(old_array, in_vecs=None, out_vecs=None, verbose=True):
        if (verbose): print('*** rotateSymmetricMatrix ***')
        
        nb_cells = old_array.GetNumberOfTuples()
        new_array = createFloatArray("", 6, nb_cells)
        
        for num_cell in range(nb_cells):
            old_matrix = vec_col_to_mat_sym(old_array.GetTuple(num_cell))
        
            if (in_vecs == None):
                in_R = numpy.eye(3)
            else:
                in_R = numpy.transpose(numpy.array([in_vecs[0].GetTuple(num_cell),
                                                    in_vecs[1].GetTuple(num_cell),
                                                    in_vecs[2].GetTuple(num_cell)]))
        
            if (out_vecs == None):
                out_R = numpy.eye(3)
            else:
                out_R = numpy.transpose(numpy.array([out_vecs[0].GetTuple(num_cell),
                                                     out_vecs[1].GetTuple(num_cell),
                                                     out_vecs[2].GetTuple(num_cell)]))
        
            R = numpy.dot(numpy.transpose(in_R), out_R)
        
            new_matrix = numpy.dot(numpy.dot(numpy.transpose(R), old_matrix), R)
        
            new_array.InsertTuple(num_cell, mat_sym_to_vec_col(new_matrix))
        
        return new_array

#def rotateSymmetricMatrixToOrFromLocalCylindricalCoordinates(mesh, old_field_name, new_field_name, to_or_from, verbose=True):
    #nb_cells = mesh.GetNumberOfCells()

    #old_field_array = mesh.GetCellData().GetArray(old_field_name)

    #eR_array = mesh.GetCellData().GetArray('eR')
    #eC_array = mesh.GetCellData().GetArray('eC')
    #eL_array = mesh.GetCellData().GetArray('eL')

    #if (verbose): print 'Rotating field and Filling mesh...'

    #new_field_array = vtk.vtkFloatArray()
    #new_field_array.SetName(new_field_name)
    #new_field_array.SetNumberOfComponents(6)
    #new_field_array.SetNumberOfTuples(nb_cells)

    #for num_cell in range(nb_cells):
        #old_matrix = vec_col_to_mat_sym(old_field_array.GetTuple(num_cell))

        #eR = eR_array.GetTuple(num_cell)
        #eC = eC_array.GetTuple(num_cell)
        #eL = eL_array.GetTuple(num_cell)

        #R = numpy.transpose(numpy.array([eR, eC, eL]))

        #if   (to_or_from ==   'to'): new_matrix = numpy.dot(numpy.dot(numpy.transpose(R), old_matrix), R)
        #elif (to_or_from == 'from'): new_matrix = numpy.dot(numpy.dot(R, old_matrix), numpy.transpose(R))

        #new_field_array.InsertTuple(num_cell, mat_sym_to_vec_col(new_matrix))

    #mesh.GetCellData().AddArray(new_field_array)

#def rotateSymmetricMatrixToLocalCylindricalCoordinates(mesh, old_field_name, new_field_name, verbose=True):
    #rotateSymmetricMatrixToOrFromLocalCylindricalCoordinates(mesh, old_field_name, new_field_name, 'to', verbose)

#def rotateSymmetricMatrixFromLocalCylindricalCoordinates(mesh, old_field_name, new_field_name, verbose=True):
    #rotateSymmetricMatrixToOrFromLocalCylindricalCoordinates(mesh, old_field_name, new_field_name, 'from', verbose)
