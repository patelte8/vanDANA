########################################################################

import sys
import numpy
import vtk
from vtk_py3 import *

########################################################################

def transform_mesh_w_4x4mat(infilename, outfilename, matrix):

        ugrid = vtk.vtkUnstructuredGrid()
        
        if(infilename[len(infilename)-4:len(infilename)] == ".vtu"):
        	ugrid = readXMLUGrid(infilename)
        elif(infilename[len(infilename)-4:len(infilename)] == ".vtk"):
        	ugrid = readUGrid(infilename)
        
        rot_mat = vtk.vtkMatrix4x4()
        rot_mat.DeepCopy(matrix)
        
        transform = vtk.vtkTransform()
        transform.SetMatrix(rot_mat)
        transform.Update()
        
        transformfilter = vtk.vtkTransformFilter()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        	transformfilter.SetInputData(ugrid)
        else:
        	transformfilter.SetInput(ugrid)
        transformfilter.SetTransform(transform)
        transformfilter.Update()
        
        if(outfilename[len(outfilename)-4:len(outfilename)] == ".vtu"):
        	writeXMLUGrid(transformfilter.GetOutput(), outfilename)
        elif(outfilename[len(outfilename)-4:len(outfilename)] == ".vtk"):
        	writeUGrid(transformfilter.GetOutput(), outfilename)
        
        
        

if (__name__ == "__main__"):
	
        print(len(sys.argv))
        assert (len(sys.argv) == 19), 'Number of arguments must be 3.'
        infilename = sys.argv[1]
        outfilename = sys.argv[2]
        matrix = map(float,sys.argv[3:19])
        
        print(matrix)
        
        transform_mesh_w_4x4mat(infilename, outfilename, matrix)	


