import vtk as vtk
from vtk_py3 import *

def extractUGridBasedOnThreshold(ugrid, arrayname, thresholdval):

        selectionNode = vtk.vtkSelectionNode();
        selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL);
        ugrid.GetCellData().SetActiveScalars(arrayname);
        
        idArray = vtk.vtkIdTypeArray()
        if(type(thresholdval) == list):
        	selectionNode.SetContentType(vtk.vtkSelectionNode.THRESHOLDS);
        	idArray.SetNumberOfTuples(len(thresholdval)*2)
        	cnt = 0
        	for tval in thresholdval:
        		idArray.SetValue(cnt, tval)
        		idArray.SetValue(cnt+1, tval)
        		cnt += 2
        else:
        	selectionNode.SetContentType(vtk.vtkSelectionNode.THRESHOLDS);
        	idArray.SetNumberOfTuples(2)
        	idArray.SetValue(0, thresholdval)
        	idArray.SetValue(1, thresholdval)
        
        selectionNode.SetSelectionList(idArray)
        
        selection = vtk.vtkSelection();
        selection.AddNode(selectionNode);
        
        extract = vtk.vtkExtractSelectedThresholds();
        extract.SetInputData(0, ugrid);
        extract.SetInputData(1, selection);
        extract.Update();
        
        
        return extract.GetOutput()


