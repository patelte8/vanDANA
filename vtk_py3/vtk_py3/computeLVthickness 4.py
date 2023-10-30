#######################################################################

import numpy
import vtk
import os
import sys
import math

import vtk_py3 as myVTK

########################################################################

def computeLVthickness(ugrid, field_type="point"):


        assert (field_type in ["point", "cell"]), "\"field_type\" must be \"point\" or \"cell\". Aborting."
        
        pdata = myVTK.convertUGridtoPdata(ugrid)
        C = myVTK.getcentroid(pdata)
        ztop = pdata.GetBounds()[5]
        C = [C[0], C[1], ztop-0.05]
        clippedheart = myVTK.clipheart(pdata, C, [0,0,1], True) 
        epi , endo= myVTK.splitDomainBetweenEndoAndEpi(clippedheart)
        
        cleanepipdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        	cleanepipdata.SetInputData(epi)
        else:
        	cleanepipdata.SetInput(epi)
        cleanepipdata.Update()
        cleanepi = cleanepipdata.GetOutput()	
        
        cleanendopdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        	cleanendopdata.SetInputData(endo)
        else:
        	cleanendopdata.SetInput(endo)
        cleanendopdata.Update()
        cleanendo = cleanendopdata.GetOutput()	
        
        if(field_type == "point"):
        
        	pointLocator = vtk.vtkPointLocator()
        	pointLocator.SetDataSet(cleanepi)
        	pointLocator.BuildLocator()
        
        	thickness = vtk.vtkFloatArray()
        	thickness.SetName("Thickness")
        	thickness.SetNumberOfComponents(1)
        
        	closest_epi_ptid = vtk.vtkIdList()
        	for ptid in range(0, cleanendo.GetNumberOfPoints()):
        		endopt = cleanendo.GetPoints().GetPoint(ptid)
        		closest_epi_ptid = pointLocator.FindClosestPoint(endopt)
        		closestepipt = cleanepi.GetPoints().GetPoint(closest_epi_ptid)
        
        		distance = math.sqrt(vtk.vtkMath().Distance2BetweenPoints(endopt, closestepipt))
        		thickness.InsertNextValue(distance)
        
        	cleanendo.GetPointData().AddArray(thickness)
        
        if(field_type == "cell"):
        	pendo_cellcenter = vtk.vtkCellCenters()
        	if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        		pendo_cellcenter.SetInputData(cleanendo)
        	else:
        		pendo_cellcenter.SetInput(cleanendo)
        	pendo_cellcenter.Update()
        	cleanendo_cellcenter = pendo_cellcenter.GetOutput()
        
        	pepi_cellcenter = vtk.vtkCellCenters()
        	if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        		pepi_cellcenter.SetInputData(cleanepi)
        	else:
        		pepi_cellcenter.SetInput(cleanepi)
        	pepi_cellcenter.Update()
        	cleanepi_cellcenter = pepi_cellcenter.GetOutput()
        
        	pointLocator = vtk.vtkPointLocator()
        	pointLocator.SetDataSet(cleanepi_cellcenter)
        	pointLocator.BuildLocator()
        
        	thickness = vtk.vtkFloatArray()
        	thickness.SetName("Thickness")
        	thickness.SetNumberOfComponents(1)
        
        	closest_epi_ptid = vtk.vtkIdList()
        	for ptid in range(0, cleanendo.GetNumberOfCells()):
        		endopt = cleanendo_cellcenter.GetPoints().GetPoint(ptid)
        		closest_epi_ptid = pointLocator.FindClosestPoint(endopt)
        		closestepipt = cleanepi_cellcenter.GetPoints().GetPoint(closest_epi_ptid)
        
        		distance = math.sqrt(vtk.vtkMath().Distance2BetweenPoints(endopt, closestepipt))
        		thickness.InsertNextValue(distance)
        
        	cleanendo.GetCellData().AddArray(thickness)
        
        
        return cleanendo

		

	


	



	
	



