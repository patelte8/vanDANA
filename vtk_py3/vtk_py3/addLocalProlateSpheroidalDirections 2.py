########################################################################

import sys
import math
import numpy
import vtk

from vtk_py3.createFloatArray               import *
from vtk_py3.getABPointsFromBoundsAndCenter import *
from vtk_py3.getCellCenters                 import *
from vtk_py3.getPDataNormals                import *
from vtk_py3.writePData			    import *

########################################################################

def addLocalProlateSpheroidalDirections(ugrid_wall,
                                        pdata_end,
                                        pdata_epi,
                                        type_of_support="cell",
					epiflip=False,
					endoflip=False,
					apexflip=False,
                                        points_AB=None,
					eCCname="eCC",
					eLLname="eLL",
					eRRname="eRR",
                                        verbose=True):

    if (verbose): print('*** addLocalProlateSpheroidalDirections ***')

    if (points_AB == None):
        points_AB = getABPointsFromBoundsAndCenter(pdata_epi, verbose)
    assert (points_AB.GetNumberOfPoints() == 2), "points_AB must have two points. Aborting."
    point_A = numpy.array([0.]*3)
    point_B = numpy.array([0.]*3)
    points_AB.GetPoint(0, point_A)
    points_AB.GetPoint(1, point_B)
    if(apexflip):
        eL  = point_A - point_B
    else:
        eL  = point_B - point_A
    eL /= numpy.linalg.norm(eL)

    if (type_of_support == "cell"):
        pdata_cell_centers = getCellCenters(ugrid_wall)

    if (verbose): print("Computing cell normals...")

    if(epiflip):
    	pdata_epi = getPDataNormals(pdata_epi, flip=1)
    else:
    	pdata_epi = getPDataNormals(pdata_epi, flip=0)

    if(endoflip):
    	pdata_end = getPDataNormals(pdata_end, flip=1)
    else:
    	pdata_end = getPDataNormals(pdata_end, flip=0)
   
    if (verbose): print("Computing surface bounds...")

    bounds_end = pdata_end.GetBounds()
    bounds_epi = pdata_epi.GetBounds()
    z_min_end = bounds_end[4]
    z_min_epi = bounds_epi[4]
    z_max_end = bounds_end[5]
    z_max_epi = bounds_epi[5]
    L_end = z_max_end-z_min_end
    L_epi = z_max_epi-z_min_epi

    if (verbose): print("Initializing cell locators...")

    cell_locator_end = vtk.vtkCellLocator()
    cell_locator_end.SetDataSet(pdata_end)
    cell_locator_end.Update()

    cell_locator_epi = vtk.vtkCellLocator()
    cell_locator_epi.SetDataSet(pdata_epi)
    cell_locator_epi.Update()

    closest_point_end = [0.]*3
    closest_point_epi = [0.]*3
    generic_cell = vtk.vtkGenericCell()
    cellId_end = vtk.mutable(0)
    cellId_epi = vtk.mutable(0)
    subId = vtk.mutable(0)
    dist_end = vtk.mutable(0.)
    dist_epi = vtk.mutable(0.)

    if (verbose): print("Computing local prolate spheroidal directions...")

    if (type_of_support == "cell"):
        nb_cells = ugrid_wall.GetNumberOfCells()
    elif (type_of_support == "point"):
        nb_cells = ugrid_wall.GetNumberOfPoints()

    farray_norm_dist_end = createFloatArray("norm_dist_end", 1, nb_cells)
    farray_norm_dist_epi = createFloatArray("norm_dist_epi", 1, nb_cells)

    farray_norm_z_end = createFloatArray("norm_z_end", 1, nb_cells)
    farray_norm_z_epi = createFloatArray("norm_z_epi", 1, nb_cells)

    farray_eRR = createFloatArray(eRRname, 3, nb_cells)
    farray_eCC = createFloatArray(eCCname, 3, nb_cells)
    farray_eLL = createFloatArray(eLLname, 3, nb_cells)

    for num_cell in range(nb_cells):
        if (type_of_support == "cell"):
            cell_center = numpy.array(pdata_cell_centers.GetPoints().GetPoint(num_cell))
        elif (type_of_support == "point"):
            cell_center = numpy.array(ugrid_wall.GetPoints().GetPoint(num_cell))
        cell_locator_end.FindClosestPoint(cell_center, closest_point_end, generic_cell, cellId_end, subId, dist_end)
        cell_locator_epi.FindClosestPoint(cell_center, closest_point_epi, generic_cell, cellId_epi, subId, dist_epi)

        norm_dist_end = dist_end/(dist_end+dist_epi)
        norm_dist_epi = dist_epi/(dist_end+dist_epi)
        farray_norm_dist_end.InsertTuple(num_cell, [norm_dist_end])
        farray_norm_dist_epi.InsertTuple(num_cell, [norm_dist_epi])

        norm_z_end = (closest_point_end[2]-z_min_end)/L_end
        norm_z_epi = (closest_point_epi[2]-z_min_epi)/L_epi
        farray_norm_z_end.InsertTuple(num_cell, [norm_z_end])
        farray_norm_z_epi.InsertTuple(num_cell, [norm_z_epi])

        normal_end = numpy.reshape(pdata_end.GetCellData().GetNormals().GetTuple(cellId_end), (3))
        normal_epi = numpy.reshape(pdata_epi.GetCellData().GetNormals().GetTuple(cellId_epi), (3))
        eRR  = -1*(1.-norm_dist_end) * normal_end + (1.-norm_dist_epi) * normal_epi
        eRR /= numpy.linalg.norm(eRR)
        eCC  = numpy.cross(eL, eRR)
        eCC /= numpy.linalg.norm(eCC)
        eLL  = numpy.cross(eRR, eCC)
        farray_eRR.InsertTuple(num_cell, eRR)
        farray_eCC.InsertTuple(num_cell, eCC)
        farray_eLL.InsertTuple(num_cell, eLL)

    if (verbose): print("Filling mesh...")

    if (type_of_support == "cell"):
        ugrid_wall.GetCellData().AddArray(farray_norm_dist_end)
        ugrid_wall.GetCellData().AddArray(farray_norm_dist_epi)
        ugrid_wall.GetCellData().AddArray(farray_norm_z_end)
        ugrid_wall.GetCellData().AddArray(farray_norm_z_epi)
        ugrid_wall.GetCellData().AddArray(farray_eRR)
        ugrid_wall.GetCellData().AddArray(farray_eCC)
        ugrid_wall.GetCellData().AddArray(farray_eLL)
    elif (type_of_support == "point"):
        ugrid_wall.GetPointData().AddArray(farray_norm_dist_end)
        ugrid_wall.GetPointData().AddArray(farray_norm_dist_epi)
        ugrid_wall.GetPointData().AddArray(farray_norm_z_end)
        ugrid_wall.GetPointData().AddArray(farray_norm_z_epi)
        ugrid_wall.GetPointData().AddArray(farray_eRR)
        ugrid_wall.GetPointData().AddArray(farray_eCC)
        ugrid_wall.GetPointData().AddArray(farray_eLL)

if (__name__ == "__main__"):
    assert (len(sys.argv) in [2]), "Number of arguments must be 1. Aborting."
    basename = sys.argv[1]
    ugrid_wall = readUGrid(basename + "-Mesh.vtk")
    pdata_end = readSTL(basename + "-End.stl")
    pdata_epi = readSTL(basename + "-Epi.stl")
    addLocalProlateSpheroidalDirections(ugrid_wall, pdata_end, pdata_epi)
    writeUGrid(ugrid_wall, basename + "-Mesh.vtk")
