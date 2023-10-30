########################################################################
import vtk
import numpy as np
import vtk_py3 as vtk_py
import sys
import inspect
import os

########################################################################

def create_BiVmesh(epicutfilename, LVendocutfilename, RVendocutfilename, casename="BiV", meshsize=0.5, gmshcmd="gmsh", iswritemesh=False, verbose=True):


        if (verbose): print('*** create_BiVmesh ***')
        
        cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
        
        geofile = cur_dir + "/BiV.geo"
        tempgeofile = "BiV_temp.geo"
        print(tempgeofile)
        print(os.getcwd())
        
        meshfilename = casename+".vtk"
        
        cmd = "cp " + geofile + " " + tempgeofile
        os.system(cmd)
        
        cmd = "sed -i.bak s/'<<Meshsize>>'/'" + str(meshsize) + "'/g " + tempgeofile
        os.system(cmd)
        cmd = "sed -i.bak s/'<<LVfilename>>'/'" + "\""+ str(LVendocutfilename).replace("/", "\/") + "\"" + "'/g " + tempgeofile
        os.system(cmd)
        cmd = "sed -i.bak s/'<<RVfilename>>'/'" + "\""+ str(RVendocutfilename).replace("/", "\/") + "\"" + "'/g " + tempgeofile
        os.system(cmd)
        cmd = "sed -i.bak s/'<<Epifilename>>'/'" + "\""+ str(epicutfilename).replace("/", "\/") + "\"" + "'/g " + tempgeofile
        os.system(cmd)
        cmd = gmshcmd+" -3 BiV_temp.geo -o " + meshfilename
        os.system(cmd)
        cmd = "rm BiV_temp.geo "
        os.system(cmd)
        
        ugrid = vtk_py.readUGrid(meshfilename)
        #cmd = "rm meshfilename"
        os.system(cmd)
        
        return ugrid





