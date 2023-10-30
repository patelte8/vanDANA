import vtk
import vtk_py3 as vtk_py
import dolfin as dolfin
import numpy as np

def extractFeNiCsBiVFacet(ugrid, geometry="BiV"):

        tol = 1e-2
        
        #ugrid = vtk_py.readUGrid(meshfilename)
        
        # Extract surface
        geom = vtk.vtkGeometryFilter()
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        	geom.SetInput(ugrid)
        else:
        	geom.SetInputData(ugrid)
        geom.Update()
        surf = geom.GetOutput()
        
        bc_pts_locator = []
        bc_pts = []
        bc_pts_range = []
        bc_pts_map = []
        
        # Extract Surface Normal
        normal = vtk.vtkPolyDataNormals()
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        	normal.SetInput(surf)
        else:
        	normal.SetInputData(surf)
        normal.ComputeCellNormalsOn()
        normal.Update()
        surf_w_norm = normal.GetOutput()
        
        #vtk_py.writePData(normal.GetOutput(), "normal.vtk")
        
        zmax = surf_w_norm.GetBounds()[5]
        
        surf_w_norm.BuildLinks()
        idlist = vtk.vtkIdList()
        basecellidlist = vtk.vtkIdTypeArray()
        basesurf = vtk.vtkPolyData()
        for p in range(0, surf_w_norm.GetNumberOfCells()):
        	zvec = surf_w_norm.GetCellData().GetNormals().GetTuple3(p)[2]
        
        	surf_w_norm.GetCellPoints(p, idlist)
        	zpos = surf_w_norm.GetPoints().GetPoint(idlist.GetId(0))[2]
        
        	if((abs(zvec - 1.0) < tol or abs(zvec + 1.0) < tol) and (abs(zmax - zpos) < tol)):
        		surf_w_norm.DeleteCell(p)
        		basecellidlist.InsertNextValue(p)
        
        basesurf = vtk_py.extractCellFromPData(basecellidlist, surf)
        baseptlocator = vtk.vtkPointLocator()
        baseptlocator.SetDataSet(basesurf)
        baseptlocator.BuildLocator()
        
        #######################################################################
        
        surf_w_norm.RemoveDeletedCells()
        
        
        cleanpdata = vtk.vtkCleanPolyData()
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        	cleanpdata.SetInput(surf_w_norm)
        else:
        	cleanpdata.SetInputData(surf_w_norm)
        cleanpdata.Update()
        
        connfilter = vtk.vtkPolyDataConnectivityFilter()
        if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        	connfilter.SetInput(cleanpdata.GetOutput())
        else:
        	connfilter.SetInputData(cleanpdata.GetOutput())
        connfilter.Update()
        
        print("Total_num_points = ",  cleanpdata.GetOutput().GetNumberOfPoints())
        tpt = 0
        
        if(geometry=="BiV"):
        	nsurf = 3
        else:
        	nsurf = 2
        
        
        for p in range(0,nsurf):
        
        	pts = vtk.vtkPolyData()
        
        	connfilter.SetExtractionModeToSpecifiedRegions()
        	[connfilter.DeleteSpecifiedRegion(k) for k in range(0,nsurf)]
        	connfilter.AddSpecifiedRegion(p)
        	connfilter.ScalarConnectivityOff()
        	connfilter.FullScalarConnectivityOff()
        	connfilter.Update()
        
        	cleanpdata2 = vtk.vtkCleanPolyData()
        	if(vtk.vtkVersion().GetVTKMajorVersion() < 6):
        		cleanpdata2.SetInput(connfilter.GetOutput())
        	else:
        		cleanpdata2.SetInputData(connfilter.GetOutput())
        	cleanpdata2.Update()
        
        	pts.DeepCopy(cleanpdata2.GetOutput())
        
        	tpt = tpt + cleanpdata2.GetOutput().GetNumberOfPoints()
        
        	ptlocator = vtk.vtkPointLocator()
        	ptlocator.SetDataSet(pts)
        	ptlocator.BuildLocator()
        
        	bc_pts_locator.append(ptlocator)
        	bc_pts.append(pts)
        	bc_pts_range.append([abs(pts.GetBounds()[k+1] - pts.GetBounds()[k]) for k in range(0, 6, 2)])
        
        
        #vtk_py.writePData(connfilter.GetOutput(), "/home/likchuan/Research/fenicsheartmesh/ellipsoidal/Geometry/test.vtk")
        
        print("Total_num_points = ",  tpt)
        
        Epiid = np.argmax(np.array([max(pts) for pts in bc_pts_range]))
        maxzrank =  np.array([pts[2] for pts in bc_pts_range]).argsort()
        
        
        if(geometry=="BiV"):
        	LVid = maxzrank[1] 
        	RVid = 3 - (LVid + Epiid)
        	bc_pts_map = [4, 4, 4, 4]
        	bc_pts_map[Epiid] = 1; bc_pts_map[LVid] = 2; bc_pts_map[RVid] = 3
        	baseid  = 3;
        else:
        	LVid = maxzrank[0]
        	bc_pts_map = [4, 4, 4]
        	bc_pts_map[Epiid] = 1; bc_pts_map[LVid] = 2
        	baseid  = 2;
        
        
        bc_pts_locator.append(baseptlocator)
        bc_pts.append(basesurf)
        
        
        dolfin_mesh = vtk_py.convertUGridToXMLMesh(ugrid)
        dolfin_facets = dolfin.MeshFunction('size_t', dolfin_mesh, 2)
        dolfin_facets.set_all(0)
        
        for facet in dolfin.SubsetIterator(dolfin_facets, 0):
        	for locator in range(0,nsurf+1):
        		cnt = 0
        		for p in range(0,3):
        			v0 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(0)
        			v1 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(1)
        			v2 =  dolfin.Vertex(dolfin_mesh, facet.entities(0)[p]).x(2)
        			ptid = bc_pts_locator[locator].FindClosestPoint(v0, v1, v2)
        			x0 =  bc_pts[locator].GetPoints().GetPoint(ptid)
        			dist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], x0)
        			if(dist < 1e-5):
        				cnt = cnt + 1
        		if(cnt == 3):
        			dolfin_facets[facet] = bc_pts_map[locator]
        				
        
        dolfin_edges = dolfin.MeshFunction('size_t', dolfin_mesh, 1)
        dolfin_edges.set_all(0)
        
        epilocator = Epiid
        lvendolocator = LVid
        
        for edge in dolfin.SubsetIterator(dolfin_edges, 0):
        	cnt_epi = 0; cnt_lvendo = 0;
        	for p in range(0,2):
        		v0 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(0)
        		v1 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(1)
        		v2 =  dolfin.Vertex(dolfin_mesh, edge.entities(0)[p]).x(2)
        
        		epiptid = bc_pts_locator[epilocator].FindClosestPoint(v0, v1, v2)
        		epix0 =  bc_pts[epilocator].GetPoints().GetPoint(epiptid)
        		epidist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], epix0)
        
        		topptid = bc_pts_locator[baseid].FindClosestPoint(v0, v1, v2)
        		topx0 =  bc_pts[baseid].GetPoints().GetPoint(topptid)
        		topdist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], topx0)
        
        		lvendoptid = bc_pts_locator[lvendolocator].FindClosestPoint(v0, v1, v2)
        		lvendox0 =  bc_pts[lvendolocator].GetPoints().GetPoint(lvendoptid)
        		lvendodist = vtk.vtkMath.Distance2BetweenPoints([v0,v1,v2], lvendox0)
        
        		if(topdist < 1e-5 and epidist < 1e-5):
        			cnt_epi = cnt_epi + 1
        
        		if(topdist < 1e-5 and lvendodist < 1e-5):
        			cnt_lvendo = cnt_lvendo + 1
        
        		if(cnt_epi == 2):
        			dolfin_edges[edge] = 1
        
        		if(cnt_lvendo == 2):
        			dolfin_edges[edge] = 2
        
        
        
        return dolfin_mesh, dolfin_facets, dolfin_edges	
