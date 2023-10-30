import vtk
import vtk_py3 as vtk_py

def extractCellFromPData(cellidlist, pdata):
        
        selectionNode = vtk.vtkSelectionNode()
        selectionNode.SetFieldType(vtk.vtkSelectionNode.CELL);
        selectionNode.SetContentType(vtk.vtkSelectionNode.INDICES);
        selectionNode.SetSelectionList(cellidlist);
        selection = vtk.vtkSelection();
        selection.AddNode(selectionNode);
        extractSelection = vtk.vtkExtractSelection();
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        	extractSelection.SetInput(0, pdata);
        	extractSelection.SetInput(1, selection);
        else:	
        	extractSelection.SetInputData(0, pdata);
        	extractSelection.SetInputData(1, selection);
        extractSelection.Update();
        extractbase = vtk.vtkGeometryFilter()
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        	extractbase.SetInput(extractSelection.GetOutput())
        else:
        	extractbase.SetInputData(extractSelection.GetOutput())
        extractbase.Update()
        
        return extractbase.GetOutput()


