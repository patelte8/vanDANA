########################################################################

import math
import numpy
import vtk

from vtk_py3.createFloatArray               import *
from vtk_py3.getABPointsFromBoundsAndCenter import *
from vtk_py3.getCellCenters                 import *

########################################################################

def addLocalCylindricalDirections(ugrid_wall,
                                  type_of_support="cell",
                                  points_AB=None,
                                  verbose=True):

        if (verbose): print("*** addLocalCylindricalDirections ***")
        
        if (points_AB == None):
            points_AB = getABPointsFromBoundsAndCenter(ugrid_wall, verbose)
            #print points_AB
        assert (points_AB.GetNumberOfPoints() >= 2), "points_AB must have at least two points. Aborting."
        point_A = numpy.array([0.]*3)
        point_B = numpy.array([0.]*3)
        points_AB.GetPoint(                              0, point_A)
        points_AB.GetPoint(points_AB.GetNumberOfPoints()-1, point_B)
        #if (verbose): print "point_A =", point_A
        #if (verbose): print "point_B =", point_B
        eL  = point_B - point_A
        eL /= numpy.linalg.norm(eL)
        #if (verbose): print "eL =", eL
        
        if (type_of_support == "cell"):
            pdata_cell_centers = getCellCenters(ugrid_wall)
        
        if (type_of_support == "cell"):
            nb_cells = ugrid_wall.GetNumberOfCells()
        elif (type_of_support == "point"):
            nb_cells = ugrid_wall.GetNumberOfPoints()
        
        farray_eR = createFloatArray("eRR", 3, nb_cells)
        farray_eC = createFloatArray("eCC", 3, nb_cells)
        farray_eL = createFloatArray("eLL", 3, nb_cells)
        
        farray_r = createFloatArray("r", 1, nb_cells)
        farray_t = createFloatArray("t", 1, nb_cells)
        farray_z = createFloatArray("z", 1, nb_cells)
        
        for num_cell in range(nb_cells):
            #if (verbose): print "num_cell =", num_cell
        
            if (type_of_support == "cell"):
                cell_center = numpy.array(pdata_cell_centers.GetPoints().GetPoint(num_cell))
            elif (type_of_support == "point"):
                cell_center = numpy.array(ugrid_wall.GetPoints().GetPoint(num_cell))
        
            #if (verbose): print "cell_center =", cell_center
        
            #eR  = cell_center - point_A
            #eR -= numpy.dot(eR,eL) * eL
            #eR /= numpy.linalg.norm(eR)
        
            #eC  = numpy.cross(eL, eR)
        
            eR  = cell_center - point_A
            eC  = numpy.cross(eL, eR)
            eC /= numpy.linalg.norm(eC)
            eR  = numpy.cross(eC, eL)
        
            if (numpy.dot(eR,eC) > 1e-6) or (numpy.dot(eR,eL) > 1e-6) or (numpy.dot(eC,eL) > 1e-6): print("WTF?!")
        
            farray_eR.InsertTuple(num_cell, eR)
            farray_eC.InsertTuple(num_cell, eC)
            farray_eL.InsertTuple(num_cell, eL)
        
            r = numpy.dot(cell_center - point_A, eR)
            farray_r.InsertTuple(num_cell, [r])
        
            t  = math.atan2(eR[1], eR[0])
            t += (t<0.)*(2*math.pi)
            farray_t.InsertTuple(num_cell, [t])
        
            z = numpy.dot(cell_center - point_A, eL)
            farray_z.InsertTuple(num_cell, [z])
        
        if (type_of_support == "cell"):
            ugrid_wall.GetCellData().AddArray(farray_eR)
            ugrid_wall.GetCellData().AddArray(farray_eC)
            ugrid_wall.GetCellData().AddArray(farray_eL)
            ugrid_wall.GetCellData().AddArray(farray_r)
            ugrid_wall.GetCellData().AddArray(farray_t)
            ugrid_wall.GetCellData().AddArray(farray_z)
        elif (type_of_support == "point"):
            ugrid_wall.GetPointData().AddArray(farray_eR)
            ugrid_wall.GetPointData().AddArray(farray_eC)
            ugrid_wall.GetPointData().AddArray(farray_eL)
            ugrid_wall.GetPointData().AddArray(farray_r)
            ugrid_wall.GetPointData().AddArray(farray_t)
            ugrid_wall.GetPointData().AddArray(farray_z)
        
        return ugrid_wall

if (__name__ == "__main__"):
    assert (len(sys.argv) in [2]), "Number of arguments must be 1. Aborting."
    writeUGrid(addLocalCylindricalDirections(readUGrid(sys.argv[1])), sys.argv[1])

