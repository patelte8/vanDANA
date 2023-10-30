########################################################################

import sys
import vtk
import math

#from mat_vec_tools import *
from vtk_py3.clean_pdata import *

########################################################################

def extract_slice_thickness_LVonly(nsubdivision, domain, verbose=True):

        refinedmesh = vtk.vtkLoopSubdivisionFilter()
        refinedmesh.SetNumberOfSubdivisions(nsubdivision);
        refinedmesh.SetInput(domain)
        refinedmesh.Update()
        
        bds = domain.GetBounds()
        ctr = [0,0,0]
        ctr[0] = 0.5*(bds[0] + bds[1])
        ctr[1] = 0.5*(bds[2] + bds[3])
        ctr[2] = 0.5*(bds[4] + bds[5])
        
        featureEdges = vtk.vtkFeatureEdges();
        featureEdges.SetInput(refinedmesh.GetOutput());
        featureEdges.BoundaryEdgesOn();
        featureEdges.FeatureEdgesOff();
        featureEdges.ManifoldEdgesOff();
        featureEdges.NonManifoldEdgesOff();
        featureEdges.Update();
        
        connectfilter = vtk.vtkPolyDataConnectivityFilter()
        connectfilter.SetInput(featureEdges.GetOutput())
        connectfilter.SetExtractionModeToSpecifiedRegions()
        connectfilter.Update()
        
        epi = vtk.vtkPolyData()
        LVendo = vtk.vtkPolyData()
        
        connectfilter.SetExtractionModeToClosestPointRegion()
        connectfilter.SetClosestPoint(ctr[0], ctr[1]+1000, ctr[2]+1000)
        connectfilter.Update()
        epi.DeepCopy(connectfilter.GetOutput())
        epi = clean_pdata(epi)
        
        
        connectfilter.SetClosestPoint(ctr[0], ctr[1], ctr[2])
        connectfilter.Update()
        LVendo.DeepCopy(connectfilter.GetOutput())
        LVendo = clean_pdata(LVendo)
        
        
        epipointlocator = vtk.vtkPointLocator()
        epipointlocator.SetDataSet(epi)
        epipointlocator.BuildLocator()
        
        LVendopointlocator = vtk.vtkPointLocator()
        LVendopointlocator.SetDataSet(LVendo)
        LVendopointlocator.BuildLocator()
        
        epidistance = vtk.vtkFloatArray()
        epidistance.SetName("Thickness")
        
        
        for p in range(0, epi.GetNumberOfPoints()):
        	pt = epi.GetPoints().GetPoint(p)
        
        	LVendoid = LVendopointlocator.FindClosestPoint(pt)
        	LVendopt = LVendo.GetPoints().GetPoint(LVendoid)
        	
        	disttoLVendo = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(pt, LVendopt))
        
        	epidistance.InsertNextValue(disttoLVendo)
        
        	
        epi.GetPointData().SetActiveScalars("Thickness")
        epi.GetPointData().SetScalars(epidistance)
        
        return epi


	


