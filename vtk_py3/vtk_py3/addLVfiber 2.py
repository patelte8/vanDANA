from dolfin import *
import vtk as vtk
from vtk_py3 import *

def addLVfiber(mesh, V, casename, endo_angle, epi_angle,  casedir, isepiflip, isendoflip, isapexflip=False):


        fiberV = Function(V)
        sheetV = Function(V)
        sheetnormV = Function(V)
        cV = Function(V)
        lV = Function(V)
        rV = Function(V)
        
        ugrid=vtk_py.convertXMLMeshToUGrid(mesh)
        pdata = vtk_py.convertUGridtoPdata(ugrid)
        C = vtk_py.getcentroid(pdata)
        if(isapexflip):
            ztop = pdata.GetBounds()[4]
            C = [C[0], C[1], ztop+0.05]
            clippedheart = vtk_py.clipheart(pdata, C, [0,0,-1], True)
        else:
            ztop = pdata.GetBounds()[5]
            C = [C[0], C[1], ztop-0.05]
            clippedheart = vtk_py.clipheart(pdata, C, [0,0,1], True)
        epi, endo= vtk_py.splitDomainBetweenEndoAndEpi(clippedheart)
        
        cleanepipdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
                cleanepipdata.SetInputData(epi)
        else:
                cleanepipdata.SetInput(epi)
        cleanepipdata.Update()
        pdata_epi = cleanepipdata.GetOutput()
        
        
        cleanendopdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
                cleanendopdata.SetInputData(endo)
        else:
                cleanendopdata.SetInput(endo)
        cleanendopdata.Update()
        pdata_endo = cleanendopdata.GetOutput()
        
        L_epi = pdata_epi.GetBounds()[5]  -  pdata_epi.GetBounds()[4]
        L_endo = pdata_endo.GetBounds()[5] - pdata_endo.GetBounds()[4]
        
        
        
        if(L_endo > L_epi):
        	pdata_epi = temp
        	pdata_epi = pdata_endo
        	pdata_endo = temp
        	
        
        
        # Quad points
        gdim = mesh.geometry().dim()
        xdofmap = V.sub(0).dofmap().dofs()
        ydofmap = V.sub(1).dofmap().dofs()
        zdofmap = V.sub(2).dofmap().dofs()
        
        if(dolfin.__version__ != '1.6.0'):
        	xq = V.tabulate_dof_coordinates().reshape((-1, gdim))
        	xq0 = xq[xdofmap]  
        else:
        	xq = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, gdim))
        	xq0 = xq[xdofmap]  
        
        # Create an unstructured grid of Gauss Points
        points = vtk.vtkPoints()
        vertices = vtk.vtkCellArray()
        ugrid = vtk.vtkUnstructuredGrid()
        cnt = 0;
        for pt in xq0:
        	points.InsertNextPoint([pt[0], pt[1], pt[2]])
        	vertex = vtk.vtkVertex()
        	vertex.GetPointIds().SetId(0, cnt)
        	vertices.InsertNextCell(vertex)
        	cnt += 1
        
        ugrid.SetPoints(points)
        ugrid.SetCells(0, vertices)
        
        CreateVertexFromPoint(ugrid)
        addLocalProlateSpheroidalDirections(ugrid, pdata_endo, pdata_epi, type_of_support="cell", epiflip=isepiflip, endoflip=isendoflip, apexflip=isapexflip)
        addLocalFiberOrientation(ugrid, endo_angle, epi_angle)
        
        fiber_vector =  ugrid.GetCellData().GetArray("fiber vectors")
        sheet_vector =  ugrid.GetCellData().GetArray("sheet vectors")
        sheetnorm_vector =  ugrid.GetCellData().GetArray("sheet normal vectors")
        
        eCC_vector =  ugrid.GetCellData().GetArray("eCC")
        eLL_vector =  ugrid.GetCellData().GetArray("eLL")
        eRR_vector =  ugrid.GetCellData().GetArray("eRR")
        
        
        
        cnt = 0
        for pt in xq0:
        
        	fvec = fiber_vector.GetTuple(cnt)
        	svec = sheet_vector.GetTuple(cnt)
        	nvec = sheetnorm_vector.GetTuple(cnt)
        
        	cvec = eCC_vector.GetTuple(cnt)
        	lvec = eLL_vector.GetTuple(cnt)
        	rvec = eRR_vector.GetTuple(cnt)
        
        	fvecnorm = sqrt(fvec[0]**2 + fvec[1]**2 + fvec[2]**2)
        	svecnorm = sqrt(svec[0]**2 + svec[1]**2 + svec[2]**2)
        	nvecnorm = sqrt(nvec[0]**2 + nvec[1]**2 + nvec[2]**2)
        
        
        	if(abs(fvecnorm - 1.0) > 1e-7 or  abs(svecnorm - 1.0) > 1e-6 or abs(nvecnorm - 1.0) > 1e-7):
        		print(fvecnorm)
        		print(svecnorm)
        		print(nvecnorm)
        
        	#print xdofmap[cnt], ydofmap[cnt], zdofmap[cnt]
        	fiberV.vector()[xdofmap[cnt]] = fvec[0]; fiberV.vector()[ydofmap[cnt]] = fvec[1]; fiberV.vector()[zdofmap[cnt]] = fvec[2];
        	sheetV.vector()[xdofmap[cnt]] = svec[0]; sheetV.vector()[ydofmap[cnt]] = svec[1]; sheetV.vector()[zdofmap[cnt]] = svec[2];
        	sheetnormV.vector()[xdofmap[cnt]] = nvec[0]; sheetnormV.vector()[ydofmap[cnt]] = nvec[1]; sheetnormV.vector()[zdofmap[cnt]] = nvec[2];
        
        	cV.vector()[xdofmap[cnt]] = cvec[0];  cV.vector()[ydofmap[cnt]] = cvec[1]; cV.vector()[zdofmap[cnt]] = cvec[2]; 
        	lV.vector()[xdofmap[cnt]] = lvec[0];  lV.vector()[ydofmap[cnt]] = lvec[1]; lV.vector()[zdofmap[cnt]] = lvec[2]; 
        	rV.vector()[xdofmap[cnt]] = rvec[0];  rV.vector()[ydofmap[cnt]] = rvec[1]; rV.vector()[zdofmap[cnt]] = rvec[2]; 
        
        
        	cnt += 1
        
        
        writeXMLUGrid(ugrid, "fiber.vtu")
        
        return fiberV, sheetV, sheetnormV, cV, lV, rV


