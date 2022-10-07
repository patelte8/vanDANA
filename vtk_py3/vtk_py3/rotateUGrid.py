########################################################################

import numpy
import vtk

########################################################################

def rotateUGrid(ugrid, rx=0.0, ry=0.0, rz=0.0, sx=1.0, sy=1.0, sz=1.0, verbose=True):
        
        transmat = vtk.vtkTransform()
        transmat.RotateX(rx)
        transmat.RotateY(ry)
        transmat.RotateZ(rz)
        transmat.Scale(sx, sy, sz)
        
        transmat.Update()
        
        tf = vtk.vtkTransformFilter()
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        	tf.SetInput(ugrid)    
        else:
        	tf.SetInputData(ugrid)    
        tf.SetTransform(transmat)
        tf.Update()
        
        return tf.GetOutput()
  

	

