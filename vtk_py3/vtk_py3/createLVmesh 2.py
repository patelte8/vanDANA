########################################################################

import sys
import vtk
import os
import inspect
import vtk_py3

########################################################################

def createLVmesh(casename, meshsize, epifilename, endofilename, verbose=True):

    if (verbose): print('*** createLVmesh ***')

    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 

    LVgeofile = cur_dir + "/LV.geo"
    LVtempgeofile = "LVtemp.geo"
    mshfilename = casename + ".msh"
    vtkfilename = casename + ".vtk"

    cmd = "cp " + LVgeofile + " " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<mesh_d>>'/'" + str(meshsize) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<Endofilename>>'/'" + endofilename + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<Epifilename>>'/'" + epifilename + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "gmsh -3 LVtemp.geo -o " + mshfilename
    os.system(cmd)
    cmd = "gmsh -3 LVtemp.geo -o " + vtkfilename
    os.system(cmd)
    cmd = "rm LVtemp.geo"
    os.system(cmd)



 
