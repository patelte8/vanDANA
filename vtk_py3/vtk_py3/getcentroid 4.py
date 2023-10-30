########################################################################

import sys
import vtk

########################################################################

def getcentroid(domain, verbose=True):

    if (verbose): print('*** Get Centroid ***')

    bds = domain.GetBounds()

    return [0.5*(bds[0]+bds[1]), 0.5*(bds[2]+bds[3]), 0.5*(bds[4]+bds[5]) ]


