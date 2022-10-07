########################################################################

import vtk

########################################################################

def getABPointsFromBoundsAndCenter(ugrid_wall,
                                   verbose=True):

        if (verbose): print('*** getABPointsFromBoundsAndCenter ***')
        
        # Note that it is assumed here that the ventricle is vertical in the global coordinates system
        bounds = ugrid_wall.GetBounds()
        center = ugrid_wall.GetCenter()
        #print "bounds =", bounds
        #print "center =", center
        point_A = [center[0], center[1], bounds[4]]
        point_B = [center[0], center[1], bounds[5]]
        #print "point_A =", point_A
        #print "point_B =", point_B
        points_AB = vtk.vtkPoints()
        points_AB.InsertNextPoint(point_A)
        points_AB.InsertNextPoint(point_B)
        #print points_AB
        return points_AB
