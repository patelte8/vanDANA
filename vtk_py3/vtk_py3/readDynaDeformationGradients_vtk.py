########################################################################

import numpy
import vtk

from createFloatArray import *

########################################################################

def readDynaDeformationGradients_vtk(mesh,
                                 hystory_files_basename,
                                 array_name, ES_timept,
                                 verbose=True):

    if (verbose): print '*** readDynaDeformationGradients ***'

    nb_cells = mesh.GetNumberOfCells()

    reader = vtk.vtkLSDynaReader()
    reader.SetDatabaseDirectory('./')
    reader.SetTimeStep(ES_timept)
    reader.Update()

    print nb_cells, reader.GetNumberOfSolidCells()

    F11 = vtk.vtkDoubleArray()
    F12 = vtk.vtkDoubleArray()
    F13 = vtk.vtkDoubleArray()
    F21 = vtk.vtkDoubleArray()
    F22 = vtk.vtkDoubleArray()
    F23 = vtk.vtkDoubleArray()
    F31 = vtk.vtkDoubleArray()
    F32 = vtk.vtkDoubleArray()
    F33 = vtk.vtkDoubleArray()



    numtuples =  reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetNumberOfTuples()
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 10, 10, F11)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 11, 11, F12)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 12, 12, F13)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 13, 13, F21)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 14, 14, F22)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 15, 15, F23)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 16, 16, F31)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 17, 17, F32)
    reader.GetOutput().GetBlock(0).GetCellData().GetArray('IntPtData').GetData(0, numtuples-1, 18, 18, F33)

    F_array = createFloatArray(array_name, 9, nb_cells)

    for num_cell in range(nb_cells):
	F_list = [F11.GetValue(num_cell), F12.GetValue(num_cell), F13.GetValue(num_cell), 
		  F21.GetValue(num_cell), F22.GetValue(num_cell), F23.GetValue(num_cell), 
		  F31.GetValue(num_cell), F32.GetValue(num_cell), F33.GetValue(num_cell)]
        F_array.InsertTuple(num_cell, F_list)


    if (verbose): print "nb_tuples = " + str(F_array.GetNumberOfTuples())

    mesh.GetCellData().AddArray(F_array)




