#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import vtk

import vtk_py3 as myVTK

########################################################################

def getCellLocator(
        mesh,
        verbose=1):

    #myVTK.myPrint(verbose, "*** getCellLocator ***")
    if (verbose): print("*** getCellLocator ***")

    cell_locator = vtk.vtkCellLocator()
    cell_locator.SetDataSet(mesh)
    cell_locator.Update()

    closest_point = [0.]*3
    generic_cell = vtk.vtkGenericCell()
    k_cell = vtk.mutable(0)
    subId = vtk.mutable(0)
    dist = vtk.mutable(0.)

    return (cell_locator,
            closest_point,
            generic_cell,
            k_cell,
            subId,
            dist)
