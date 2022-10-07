#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import math
import numpy
import vtk
import os
import sys

import vtk_py3 as myVTK

########################################################################

def computeRegionsForBiV(
        points,
        pdata_endLV,
        pdata_endRV,
        pdata_epi,
        verbose=0):

    if (verbose): print("*** computeRegionsForBiV ***")
    if (verbose): print("Initializing cell locators...")


    (cell_locator_endLV,
     closest_point_endLV,
     generic_cell,
     cellId_endLV,
     subId,
     dist_endLV) = myVTK.getCellLocator(
         mesh=pdata_endLV,
         verbose=verbose-1)
    (cell_locator_endRV,
     closest_point_endRV,
     generic_cell,
     cellId_endRV,
     subId,
     dist_endRV) = myVTK.getCellLocator(
         mesh=pdata_endRV,
         verbose=verbose-1)
    (cell_locator_epi,
     closest_point_epi,
     generic_cell,
     cellId_epi,
     subId,
     dist_epi) = myVTK.getCellLocator(
         mesh=pdata_epi,
         verbose=verbose-1)

    n_points = points.GetNumberOfPoints()

    iarray_region = myVTK.createIntArray("region_id", 1, n_points)

    for k_point in range(n_points):
        point = numpy.array(points.GetPoint(k_point))
        cell_locator_endLV.FindClosestPoint(
            point,
            closest_point_endLV,
            generic_cell,
            cellId_endLV,
            subId,
            dist_endLV)
        cell_locator_endRV.FindClosestPoint(
            point,
            closest_point_endRV,
            generic_cell,
            cellId_endRV,
            subId,
            dist_endRV)
        cell_locator_epi.FindClosestPoint(
            point,
            closest_point_epi,
            generic_cell,
            cellId_epi,
            subId,
            dist_epi)

        if   (dist_endRV == max(dist_endLV, dist_endRV, dist_epi)):
            iarray_region.SetTuple(k_point, [0])
        elif (dist_epi == max(dist_endLV, dist_endRV, dist_epi)):
            iarray_region.SetTuple(k_point, [1])
        elif (dist_endLV == max(dist_endLV, dist_endRV, dist_epi)):
            iarray_region.SetTuple(k_point, [2])

    return iarray_region

########################################################################

def addRegionsToBiV(
        ugrid_mesh,
        pdata_endLV,
        pdata_endRV,
        pdata_epi,
        verbose=0):

    #myVTK.myPrint(verbose, "*** addRegionsToBiV ***")
    if (verbose): print("*** addRegionsToBiV ***")

    points = ugrid_mesh.GetPoints()
    iarray_region = computeRegionsForBiV(
        points=points,
        pdata_endLV=pdata_endLV,
        pdata_endRV=pdata_endRV,
        pdata_epi=pdata_epi,
        verbose=verbose-1)
    ugrid_mesh.GetPointData().AddArray(iarray_region)

    cell_centers = myVTK.getCellCenters(
        mesh=ugrid_mesh,
        verbose=verbose-1)
    iarray_region = computeRegionsForBiV(
        points=cell_centers,
        pdata_endLV=pdata_endLV,
        pdata_endRV=pdata_endRV,
        pdata_epi=pdata_epi,
        verbose=verbose-1)
    ugrid_mesh.GetCellData().AddArray(iarray_region)

########################################################################

def addRegionsToBiV2D(
        ugrid_mesh,
        LVplane,
        LVpoint,
        RVplane,
	RVpoint,
        verbose=0):

    #myVTK.myPrint(verbose, "*** addRegionsToBiV2D ***")
    if (verbose): print("*** addRegionsToBiV2D ***")

    tol = 1e-5
 
    n_cells = ugrid_mesh.GetNumberOfCells()

    meshCellCenter = myVTK.getCellCenters(ugrid_mesh)
    matid = myVTK.createIntArray("region_id", 1, n_cells)

    LVplanenorm = numpy.linalg.norm(LVplane)
    LVPlane = 1.0/LVplanenorm*numpy.array(LVplane)

    RVplanenorm = numpy.linalg.norm(RVplane)
    RVplane = 1.0/RVplanenorm*numpy.array(RVplane)
		

    for ptid in range(n_cells):
        x = meshCellCenter.GetPoints().GetPoint(ptid)
        if((x[0] - LVpoint[0])*LVplane[0] + (x[1] - LVpoint[1])*LVplane[1]  > tol):
                matid.SetTuple(ptid,[0])
        elif((x[0] - RVpoint[0])*RVplane[0] + (x[1] - RVpoint[1])*RVplane[1]  > tol):
                matid.SetTuple(ptid,[1])
        else:
                matid.SetTuple(ptid,[2])

    ugrid_mesh.GetCellData().AddArray(matid)

    return ugrid_mesh

		
