########################################################################

import sys
import math
import numpy
import vtk

from vtk_py3.addLocalFiberOrientation            import *
from vtk_py3.addLocalProlateSpheroidalDirections import *
from vtk_py3.createFloatArray                    import *
from vtk_py3.getABPointsFromBoundsAndCenter      import *
from vtk_py3.readSTL                             import *
from vtk_py3.readUGrid                           import *
from vtk_py3.writeUGrid                          import *

########################################################################

def addLocalFiberOrientation_infarct(ugrid_wall,
                             		fiber_angle_end,
                             		fiber_angle_epi,
			 		inf_fiber_angle_end,
                             		inf_fiber_angle_epi,
					matid,
                             		points_AB=None,
                             		verbose=True):

    if (verbose): print('*** addLocalFiberOrientation with infarct ***')

    if (points_AB == None):
        points_AB = getABPointsFromBoundsAndCenter(ugrid_wall, verbose)
    assert (points_AB.GetNumberOfPoints() >= 2), "\"points_AB\" must have at least two points. Aborting."
    point_A = numpy.array([0.]*3)
    point_B = numpy.array([0.]*3)
    points_AB.GetPoint(                              0, point_A)
    points_AB.GetPoint(points_AB.GetNumberOfPoints()-1, point_B)
    eL  = point_B - point_A
    eL /= numpy.linalg.norm(eL)

    if (verbose): print("Computing local fiber orientation...")

    farray_norm_dist_end = ugrid_wall.GetCellData().GetArray("norm_dist_end")
    farray_norm_dist_epi = ugrid_wall.GetCellData().GetArray("norm_dist_epi")
    farray_eRR = ugrid_wall.GetCellData().GetArray("eRR")
    farray_eCC = ugrid_wall.GetCellData().GetArray("eCC")
    farray_eLL = ugrid_wall.GetCellData().GetArray("eLL")

    nb_cells = ugrid_wall.GetNumberOfCells()

    farray_fiber_angle = createFloatArray("fiber_angle", 1, nb_cells)

    farray_eF = createFloatArray("fiber vectors", 3, nb_cells)
    farray_eS = createFloatArray("sheet vectors", 3, nb_cells)
    farray_eN = createFloatArray("sheet normal vectors", 3, nb_cells)
    matid_data = createFloatArray("matid", 1, nb_cells)

    for num_cell in range(nb_cells):
        norm_dist_end = farray_norm_dist_end.GetTuple(num_cell)[0]
        norm_dist_epi = farray_norm_dist_epi.GetTuple(num_cell)[0]
        
        matid_data.InsertTuple(num_cell, [matid[num_cell]])
        if(matid[num_cell] == 1):
        	fiber_angle_in_degrees = (1.-norm_dist_end) * fiber_angle_end + (1.-norm_dist_epi) * fiber_angle_epi
        else:
        	fiber_angle_in_degrees = (1.-norm_dist_end) * inf_fiber_angle_end + (1.-norm_dist_epi) * inf_fiber_angle_epi
        
        farray_fiber_angle.InsertTuple(num_cell, [fiber_angle_in_degrees])
        
        eRR = numpy.array(farray_eRR.GetTuple(num_cell))
        eCC = numpy.array(farray_eCC.GetTuple(num_cell))
        eLL = numpy.array(farray_eLL.GetTuple(num_cell))
        
        fiber_angle_in_radians = math.pi*fiber_angle_in_degrees/180
        eF = math.cos(fiber_angle_in_radians) * eCC + math.sin(fiber_angle_in_radians) * eLL
        eS = eRR
        eN = numpy.cross(eF, eS)
        farray_eF.InsertTuple(num_cell, eF)
        farray_eS.InsertTuple(num_cell, eS)
        farray_eN.InsertTuple(num_cell, eN)

    if (verbose): print("Filling mesh...")

    ugrid_wall.GetCellData().AddArray(farray_fiber_angle)
    ugrid_wall.GetCellData().AddArray(farray_eF)
    ugrid_wall.GetCellData().AddArray(farray_eS)
    ugrid_wall.GetCellData().AddArray(farray_eN)
    ugrid_wall.GetCellData().AddArray(matid_data)

if (__name__ == "__main__"):
    assert (len(sys.argv) in [4]), "Number of arguments must be 3. Aborting."
    basename = sys.argv[1]
    ugrid_wall = readUGrid(basename + "-Mesh.vtk")
    pdata_end = readSTL(basename + "-End.stl")
    pdata_epi = readSTL(basename + "-Epi.stl")
    angle_end = float(sys.argv[2])
    angle_epi = float(sys.argv[3])
    addLocalProlateSpheroidalDirections(ugrid_wall, pdata_end, pdata_epi)
    addLocalFiberOrientation(ugrid_wall, angle_end, angle_epi)
    writeUGrid(ugrid_wall, basename + "-Mesh.vtk")
