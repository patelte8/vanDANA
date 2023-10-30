########################################################################

import sys
import math
import random
import numpy
import vtk

from createFloatArray               import *
from getABPointsFromBoundsAndCenter import *

########################################################################

def addLocalFiberOrientation2(ugrid_wall,
                              angles_end,
                              angles_epi,
                              points_AB=None,
                              sigma=0.,
                              verbose=True):

    if (verbose): print '*** addLocalFiberOrientation2 ***'

    assert (len(angles_end) == len(angles_epi)), "angles_end and angle_epi must have same length (nb_long_nodes). Aborting."
    nb_long_nodes = len(angles_end)
    dz = 1./(nb_long_nodes-1)
    nb_circ_nodes = len(angles_end[0])
    for angles in angles_end+angles_epi:
        assert (len(angles) == nb_circ_nodes), "angles lists must have same length (nb_circ_nodes). Aborting."
    dt = 2*math.pi/nb_circ_nodes

    if (points_AB == None):
        points_AB = getABPointsFromBoundsAndCenter(ugrid_wall, verbose)
        #print points_AB
    assert (points_AB.GetNumberOfPoints() >= 2), "\"points_AB\" must have at least two points. Aborting."
    point_A = numpy.array([0.]*3)
    point_B = numpy.array([0.]*3)
    points_AB.GetPoint(                              0, point_A)
    points_AB.GetPoint(points_AB.GetNumberOfPoints()-1, point_B)
    #if (verbose): print "point_A =", point_A
    #if (verbose): print "point_B =", point_B

    eL  = point_B - point_A
    eL /= numpy.linalg.norm(eL)
    #if (verbose): print "eL =", eL

    if (verbose): print "Computing local fiber orientation..."

    farray_r = ugrid_wall.GetCellData().GetArray("r")
    farray_t = ugrid_wall.GetCellData().GetArray("t")
    farray_z = ugrid_wall.GetCellData().GetArray("z")

    farray_norm_dist_end = ugrid_wall.GetCellData().GetArray("norm_dist_end")
    farray_norm_dist_epi = ugrid_wall.GetCellData().GetArray("norm_dist_epi")

    farray_norm_z_end = ugrid_wall.GetCellData().GetArray("norm_z_end")
    farray_norm_z_epi = ugrid_wall.GetCellData().GetArray("norm_z_epi")

    farray_eRR = ugrid_wall.GetCellData().GetArray("eRR")
    farray_eCC = ugrid_wall.GetCellData().GetArray("eCC")
    farray_eLL = ugrid_wall.GetCellData().GetArray("eLL")

    nb_cells = ugrid_wall.GetNumberOfCells()

    farray_fiber_angle = createFloatArray("fiber_angle", 1, nb_cells)

    farray_eF = createFloatArray("eF", 3, nb_cells)
    farray_eS = createFloatArray("eS", 3, nb_cells)
    farray_eN = createFloatArray("eN", 3, nb_cells)

    for num_cell in range(nb_cells):
        #print "num_cell = " + str(num_cell)

        t = farray_t.GetTuple(num_cell)[0]
        i_t = int(t/dt/1.000001)
        #print "i_t = " + str(i_t)

        zeta = (t - i_t*dt) / dt
        #print "zeta = " + str(zeta)

        norm_z_end = farray_norm_z_end.GetTuple(num_cell)[0]
        norm_z_epi = farray_norm_z_epi.GetTuple(num_cell)[0]
        i_z_end = int(norm_z_end/dz/1.000001)
        i_z_epi = int(norm_z_epi/dz/1.000001)
        #print "i_z_end = " + str(i_z_end)
        #print "i_z_epi = " + str(i_z_epi)

        eta_end = (norm_z_end - i_z_end*dz) / dz
        eta_epi = (norm_z_epi - i_z_epi*dz) / dz
        #print "eta_end = " + str(eta_end)
        #print "eta_epi = " + str(eta_epi)

        t_ii_end = angles_end[i_z_end][i_t%nb_circ_nodes]
        t_ji_end = angles_end[i_z_end][(i_t+1)%nb_circ_nodes]
        t_ij_end = angles_end[(i_z_end+1)][i_t%nb_circ_nodes]
        t_jj_end = angles_end[(i_z_end+1)][(i_t+1)%nb_circ_nodes]
        t_ii_epi = angles_epi[i_z_epi][i_t%nb_circ_nodes]
        t_ji_epi = angles_epi[i_z_epi][(i_t+1)%nb_circ_nodes]
        t_ij_epi = angles_epi[(i_z_epi+1)][i_t%nb_circ_nodes]
        t_jj_epi = angles_epi[(i_z_epi+1)][(i_t+1)%nb_circ_nodes]
        #print "t_ii_end = " + str(t_ii_end)
        #print "t_ji_end = " + str(t_ji_end)
        #print "t_ij_end = " + str(t_ij_end)
        #print "t_jj_end = " + str(t_jj_end)
        #print "t_ii_epi = " + str(t_ii_epi)
        #print "t_ji_epi = " + str(t_ji_epi)
        #print "t_ij_epi = " + str(t_ij_epi)
        #print "t_jj_epi = " + str(t_jj_epi)

        fiber_angle_end = t_ii_end * (1 - zeta - eta_end + zeta*eta_end) \
                        + t_ji_end * (zeta - zeta*eta_end) \
                        + t_ij_end * (eta_end - zeta*eta_end) \
                        + t_jj_end * (zeta*eta_end)
        fiber_angle_epi = t_ii_epi * (1 - zeta - eta_epi + zeta*eta_epi) \
                        + t_ji_epi * (zeta - zeta*eta_epi) \
                        + t_ij_epi * (eta_epi - zeta*eta_epi) \
                        + t_jj_epi * (zeta*eta_epi)

        norm_dist_end = farray_norm_dist_end.GetTuple(num_cell)[0]
        norm_dist_epi = farray_norm_dist_epi.GetTuple(num_cell)[0]
        fiber_angle_in_degrees = (1.-norm_dist_end) * fiber_angle_end + (1.-norm_dist_epi) * fiber_angle_epi
        if (sigma > 0.): fiber_angle_in_degrees *= random.normalvariate(1., sigma)
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

    if (verbose): print "Filling mesh..."

    ugrid_wall.GetCellData().AddArray(farray_fiber_angle)
    ugrid_wall.GetCellData().AddArray(farray_eF)
    ugrid_wall.GetCellData().AddArray(farray_eS)
    ugrid_wall.GetCellData().AddArray(farray_eN)
