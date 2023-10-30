########################################################################

import sys
import vtk
import math

#from mat_vec_tools import *
from vtk_py3.clean_pdata import *

########################################################################

def extract_slice_thickness(nsubdivision, domain, LVctr, RVctr, Epictr, verbose=True):

        refinedmesh = vtk.vtkLoopSubdivisionFilter()
        refinedmesh.SetNumberOfSubdivisions(nsubdivision);
        refinedmesh.SetInput(domain)
        refinedmesh.Update()
        
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
        RVendo = vtk.vtkPolyData()
        
        #connectfilter.DeleteSpecifiedRegion(0)
        #connectfilter.DeleteSpecifiedRegion(1)
        #connectfilter.DeleteSpecifiedRegion(2)
        #connectfilter.AddSpecifiedRegion(0)
        connectfilter.SetExtractionModeToClosestPointRegion()
        connectfilter.SetClosestPoint(Epictr[0], Epictr[1], Epictr[2])
        connectfilter.Update()
        epi.DeepCopy(connectfilter.GetOutput())
        epi = clean_pdata(epi)
        
        #connectfilter.SetExtractionModeToSpecifiedRegions()
        #connectfilter.DeleteSpecifiedRegion(0)
        #connectfilter.AddSpecifiedRegion(1)
        connectfilter.SetClosestPoint(RVctr[0], RVctr[1], RVctr[2])
        connectfilter.Update()
        RVendo.DeepCopy(connectfilter.GetOutput())
        RVendo = clean_pdata(RVendo)
        
        #connectfilter.DeleteSpecifiedRegion(1)
        #connectfilter.AddSpecifiedRegion(2)
        connectfilter.SetClosestPoint(LVctr[0], LVctr[1], LVctr[2])
        connectfilter.Update()
        LVendo.DeepCopy(connectfilter.GetOutput())
        LVendo = clean_pdata(LVendo)
        
        
        epipointlocator = vtk.vtkPointLocator()
        epipointlocator.SetDataSet(epi)
        epipointlocator.BuildLocator()
        
        LVendopointlocator = vtk.vtkPointLocator()
        LVendopointlocator.SetDataSet(LVendo)
        LVendopointlocator.BuildLocator()
        
        RVendopointlocator = vtk.vtkPointLocator()
        RVendopointlocator.SetDataSet(RVendo)
        RVendopointlocator.BuildLocator()
        
        epidistance = vtk.vtkFloatArray()
        epidistance.SetName("Thickness")
        
        
        for p in range(0, epi.GetNumberOfPoints()):
        	pt = epi.GetPoints().GetPoint(p)
        
        	RVendoid = RVendopointlocator.FindClosestPoint(pt)
        	RVendopt = RVendo.GetPoints().GetPoint(RVendoid)
        
        	LVendoid = LVendopointlocator.FindClosestPoint(pt)
        	LVendopt = LVendo.GetPoints().GetPoint(LVendoid)
        	
        	disttoRVendo = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(pt, RVendopt))
        	disttoLVendo = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(pt, LVendopt))
        
        	epidistance.InsertNextValue(min(disttoRVendo, disttoLVendo))
        
        	
        epi.GetPointData().SetActiveScalars("Thickness")
        epi.GetPointData().SetScalars(epidistance)
        
        return epi


	


