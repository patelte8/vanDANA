########################################################################

import sys
import vtk

#from mat_vec_tools import *
from vtk_py3.readSTL       import *
from vtk_py3.writeSTL      import *

########################################################################

def clipSurfacesForFullLVMesh(endo, epi, verbose=True):
    
        if (verbose): print('*** clipSurfacesForFullLVMesh ***')
        
        endo_implicit_distance = vtk.vtkImplicitPolyDataDistance()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            endo_implicit_distance.SetInputData(endo)
        else:
            endo_implicit_distance.SetInput(endo)
        
        epi_implicit_distance = vtk.vtkImplicitPolyDataDistance()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            epi_implicit_distance.SetInputData(epi)
        else:
            epi_implicit_distance.SetInput(epi)
        
        epi_clip = vtk.vtkClipPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            epi_clip.SetInputData(epi)
        else:
            epi_clip.SetInput(epi)
        epi_clip.SetClipFunction(endo_implicit_distance)
        epi_clip.GenerateClippedOutputOn()
        epi_clip.Update()
        clipped_epi = epi_clip.GetOutput(0)
        clipped_valve = epi_clip.GetOutput(1)
        
        endo_clip = vtk.vtkClipPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            endo_clip.SetInputData(endo)
        else:
            endo_clip.SetInput(endo)
        endo_clip.SetClipFunction(epi_implicit_distance)
        endo_clip.InsideOutOn()
        endo_clip.Update()
        clipped_endo = endo_clip.GetOutput(0)
        
        return clipped_endo, clipped_epi, clipped_valve

if (__name__ == "__main__"):
        assert (len(sys.argv) in [2,3]), 'Number of arguments must be 2 or 3.'
        if (len(sys.argv) == 2):
            endo_filename = sys.argv[1] + '-EndoLV.stl'
            epi_filename = sys.argv[1] + '-EpiLV.stl'
            clipped_endo_filename = sys.argv[1] + '_FullLV-Endo.stl'
            clipped_epi_filename = sys.argv[1] + '_FullLV-Epi.stl'
            clipped_valve_filename = sys.argv[1] + '_FullLV-Valve.stl'
        elif (len(sys.argv) == 3):
            endo_filename = sys.argv[1]
            epi_filename = sys.argv[2]
            clipped_endo_filename = 'clipped_endo.stl'
            clipped_epi_filename = 'clipped_epi.stl'
            clipped_valve_filename = 'clipped_valve.stl'
        endo = readSTL(endo_filename)
        epi = readSTL(epi_filename)
        clipped_endo, clipped_epi, clipped_valve = clipSurfacesForFullLVMesh(endo, epi)
        writeSTL(clipped_endo, clipped_endo_filename)
        writeSTL(clipped_epi, clipped_epi_filename)
        writeSTL(clipped_valve, clipped_valve_filename)

