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

def getCellCenters(
        mesh,
        verbose=1):

    #myVTK.myPrint(verbose, "*** getCellCenters ***")
    print("*** getCellCenters ***")
    if (verbose): print("*** getCellCenters ***")

    filter_cell_centers = vtk.vtkCellCenters()
    if(vtk.vtkVersion().GetVTKMajorVersion() > 5):
    	filter_cell_centers.SetInputData(mesh)
    else:
    	filter_cell_centers.SetInput(mesh)
    filter_cell_centers.Update()

    return filter_cell_centers.GetOutput()
