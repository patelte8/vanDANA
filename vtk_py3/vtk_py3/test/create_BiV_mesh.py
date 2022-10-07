import dolfin as dolfin
import vtk_py3 as vtk_py
import vtk

'''
BiVEpi = vtk_py.convertUGridtoPdata(vtk_py.readUGrid("BiVEpi.vtk"))
BiVLVEndo = vtk_py.convertUGridtoPdata(vtk_py.readUGrid("BiVLVEndo.vtk"))
BiVRVEndo = vtk_py.convertUGridtoPdata(vtk_py.readUGrid("BiVRVEndo.vtk"))

# Write STL
vtk_py.writeSTL(BiVEpi, "BiVEpi.stl")
vtk_py.writeSTL(BiVLVEndo, "BiVLVEndo.stl")
vtk_py.writeSTL(BiVRVEndo, "BiVRVEndo.stl")

# Create BiVMesh
BiVUgrid = vtk_py.create_BiVmesh("BiVEpi.stl", "BiVLVEndo.stl", "BiVRVEndo.stl", "BiV", 0.5)

# Compute Regions for BiV
vtk_py.addRegionsToBiV(BiVUgrid, BiVLVEndo,  BiVRVEndo,  BiVEpi)
vtk_py.writeXMLUGrid(BiVUgrid, "BiVUGrid.vtu")
'''
ugrid = vtk_py.readXMLUGrid("BiVUGrid.vtu")

# Extract UGridBasedOnThreshold
extractugrid = vtk_py.extractUGridBasedOnThreshold(ugrid, "region_id", 1)
vtk_py.writeXMLUGrid(extractugrid, "extractBiVUGrid.vtu")

# CreateVertexFromPoint
newgrid = vtk.vtkUnstructuredGrid()
newgrid.SetPoints(ugrid.GetPoints())
extractpoints = vtk_py.CreateVertexFromPoint(newgrid)
vtk_py.writeXMLUGrid(newgrid, "extractBiVPoints.vtu")

# Extract BiV Facet
mesh, facetboundaries, edges = vtk_py.extractFeNiCsBiVFacet(ugrid, geometry="BiV")

dolfin.File("facetboundaries.pvd") << facetboundaries


# Set BiV fiber
fiber_angle_param = {"mesh": mesh,\
        "facetboundaries": facetboundaries,\
        "LV_fiber_angle": [60,-60], \
        "LV_sheet_angle": [0.1, -0.1], \
        "Septum_fiber_angle": [60,-60],\
        "Septum_sheet_angle": [0.1, -0.1],\
        "RV_fiber_angle": [60,-60],\
        "RV_sheet_angle": [0.1, -0.1],\
        "LV_matid": 0,\
        "Septum_matid": 1,\
        "RV_matid": 2,\
        #"matid": matid, #meshData["matid"],\
        "isrotatept": False,\
        "isreturn": True,\
        "outfilename": "BiV",\
        "outdirectory":"./",\
        "epiid": 1,\
        "rvid": 3,\
        "lvid": 2,\
        "degree": 4}

ef, es, en = vtk_py.SetBiVFiber_Quad_PyQ(fiber_angle_param)

