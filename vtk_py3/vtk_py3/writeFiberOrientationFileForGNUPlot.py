########################################################################

import sys
import math
import random
import numpy
import vtk

########################################################################

def writeFiberOrientationFileForGNUPlot(angles_end,
                                        angles_epi,
                                        fiber_orientation_file_name,
                                        verbose=True):

    if (verbose): print '*** writeFiberOrientationFileForGNUPlot ***'

    assert (len(angles_end) == len(angles_epi)), "angles_end and angle_epi must have same length (nb_long_nodes). Aborting."
    nb_long_nodes = len(angles_end)
    dz = 1./(nb_long_nodes-1)
    nb_circ_nodes = len(angles_end[0])
    for angles in angles_end+angles_epi:
        assert (len(angles) == nb_circ_nodes), "angles lists must have same length (nb_circ_nodes). Aborting."
    dt = 360./nb_circ_nodes
    
    fiber_orientation_file = open(fiber_orientation_file_name, 'w')
    fiber_orientation_file.write('# t ang_end ang_epi z\n')
    
    nb_t = 12
    for num_t in range(nb_t+1):
        t    = float(num_t) / nb_t * 360
        i_t  = int(t/dt/1.000001)
        zeta = (t - i_t*dt) / dt

        nb_z = 10
        for num_z in range(nb_z+1):
            z   = float(num_z) / nb_z
            i_z = int(z/dz/1.000001)
            eta = (z - i_z*dz) / dz

            t_ii_end = angles_end[i_z][i_t%nb_circ_nodes]
            t_ji_end = angles_end[i_z][(i_t+1)%nb_circ_nodes]
            t_ij_end = angles_end[(i_z+1)][i_t%nb_circ_nodes]
            t_jj_end = angles_end[(i_z+1)][(i_t+1)%nb_circ_nodes]
            t_ii_epi = angles_epi[i_z][i_t%nb_circ_nodes]
            t_ji_epi = angles_epi[i_z][(i_t+1)%nb_circ_nodes]
            t_ij_epi = angles_epi[(i_z+1)][i_t%nb_circ_nodes]
            t_jj_epi = angles_epi[(i_z+1)][(i_t+1)%nb_circ_nodes]

            fiber_angle_end = t_ii_end * (1 - zeta - eta + zeta*eta) \
                            + t_ji_end * (zeta - zeta*eta) \
                            + t_ij_end * (eta - zeta*eta) \
                            + t_jj_end * (zeta*eta)
            fiber_angle_epi = t_ii_epi * (1 - zeta - eta + zeta*eta) \
                            + t_ji_epi * (zeta - zeta*eta) \
                            + t_ij_epi * (eta - zeta*eta) \
                            + t_jj_epi * (zeta*eta)

            fiber_orientation_file.write(" ".join([str(x) for x in [t, fiber_angle_end, fiber_angle_epi, z]]) + "\n")
        
        fiber_orientation_file.write("\n")

    fiber_orientation_file.close()
