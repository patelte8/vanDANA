# Partitioning Biventricular mesh for visualizing EP measurements
# This script partition the mesh into polar segments:
# Input: ugrid = vtk file containing BiVMesh
#      : center = origin of the polar coordinate system (if none, mesh centroid will be used)
#      : xaxis = direction corresponding to 0 degree in the polar coord system
#      : nsectors = number of circular partition
#      : nz = number of longitudinal partition
#      : meas = measurements of size nz*nsectors


import vtk_py as vtk_py
import vtk
import numpy as np
import math as math

def partitionmeshforEP(ugrid, nsector=6, nz=2, meas=None, center=None, xaxis=[1,0]):

    if(center == None):
    	midx = 0.5*(ugrid.GetBounds()[0] + ugrid.GetBounds()[1])
    	midy = 0.5*(ugrid.GetBounds()[2] + ugrid.GetBounds()[3])
    	center = [midx, midy]

    
    cellcenter = vtk.vtkCellCenters()
    cellcenter.SetInputData(ugrid)
    cellcenter.Update()
    
    zpartition = np.linspace(ugrid.GetBounds()[4], ugrid.GetBounds()[5], nz+1)
    apartition = np.linspace(-math.pi, math.pi, nsector+1)
    
    regid = vtk.vtkIntArray()
    data = vtk.vtkFloatArray()

    for cellid in range(0, ugrid.GetNumberOfCells()):
     	x = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[0]
     	y = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[1]
     	z = cellcenter.GetOutput().GetPoints().GetPoint(cellid)[2]
    	
    	# Determine position in z direction
    	zloc = np.argmax(zpartition>z)
    
    	# Determine position in theta direction
    	norm = np.linalg.norm([(x - midx), (y - midy)])
    	angle = np.arctan2((y - midy)/norm, (x - midx)/norm)
    	sloc = np.argmax(apartition>angle)
    
    	regloc = (zloc-1)*nsector + sloc 
    
    	regid.InsertNextValue(regloc)
	data.InsertNextValue(meas[regloc-1])
    
    regid.SetName("Regionid")
    data.SetName("EP measurements")
    
    ugrid.GetCellData().AddArray(regid)
    ugrid.GetCellData().AddArray(data)
    
    vtk_py.writeXMLUGrid(ugrid, "/home/lclee/Downloads/test.vtu")

	

		




