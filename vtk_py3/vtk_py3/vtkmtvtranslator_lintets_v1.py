########################################################################

import argparse
import glob
from numpy import *
from sets import Set
from vtk import *
import os

import vtk_py as vtk_py

########################################################################




def extract_tetra(mesh):

	idlist = vtk.vtkIdList()
	celltype =  mesh.GetCellTypesArray()
	for p in range(0, celltype.GetNumberOfTuples()):
		if(float(celltype.GetTuple(p)[0]) == 10):
			idlist.InsertNextId(p)
	
	extracted = vtk.vtkExtractCells()
	extracted.SetCellList(idlist)
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		extracted.SetInput(mesh)
	else:
		extracted.SetInputData(mesh)
	extracted.Update()

	return extracted.GetOutput()



def writepdata(filename, pdata):

	pdatawriter = vtk.vtkXMLPolyDataWriter()
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		pdatawriter.SetInput(pdata)
	else:
		pdatawriter.SetInputData(pdata)
	pdatawriter.SetFileName(filename)
	pdatawriter.Write()



def savepoints(ugrid, nodeids, ptfilename):

	# View points
	viewpts = vtk.vtkPoints()
	viewptspdata =vtk.vtkPolyData()
	for p in range(0, nodeids.GetNumberOfIds()):
		pt = [0,0,0]
		pt = ugrid.GetPoints().GetPoint(nodeids.GetId(p))
		viewpts.InsertNextPoint(pt)

	viewptspdata.SetPoints(viewpts)

	maskPoints = vtk.vtkMaskPoints()
  	maskPoints.SetOnRatio(1); 
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
  		maskPoints.SetInput(viewptspdata);
	else:
  		maskPoints.SetInputData(viewptspdata);
	maskPoints.GenerateVerticesOn();
  	maskPoints.Update();

	writepdata(ptfilename, maskPoints.GetOutput())	
		

def extract_biven_surf_nodes(ugrid, basalnormal, tol):

	Basalsurfnodes = vtk.vtkIdList()
	Episurfnodes = vtk.vtkIdList()
	RVEndosurfnodes = vtk.vtkIdList()
	LVEndosurfnodes = vtk.vtkIdList()
	Basalepisurfnodes = vtk.vtkIdList()

	Episurf = vtk.vtkPolyData()
	LVEndosurf = vtk.vtkPolyData()
	RVEndosurf = vtk.vtkPolyData()

	Idfilter = vtk.vtkIdFilter()
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		Idfilter.SetInput(ugrid)
	else:
		Idfilter.SetInputData(ugrid)
	Idfilter.Update()

	geomfilter = vtk.vtkGeometryFilter()
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		geomfilter.SetInput(Idfilter.GetOutput())
	else:
		geomfilter.SetInputData(Idfilter.GetOutput())
	geomfilter.Update()

	cleanpdata = vtk.vtkCleanPolyData()
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		cleanpdata.SetInput(geomfilter.GetOutput())
	else:
		cleanpdata.SetInputData(geomfilter.GetOutput())
	cleanpdata.Update()

	pdatanormal = vtk.vtkPolyDataNormals()
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		pdatanormal.SetInput(cleanpdata.GetOutput())
	else:
		pdatanormal.SetInputData(cleanpdata.GetOutput())
	pdatanormal.ComputeCellNormalsOn()
	pdatanormal.Update()

	#vtk_py.writeXMLPData(pdatanormal.GetOutput(), "test.vtp")
	
	reducedsurfacemesh = vtk.vtkPolyData()
	reducedsurfacemesh.DeepCopy(pdatanormal.GetOutput())
	reducedsurfacemesh.BuildLinks()
	
	bds = reducedsurfacemesh.GetBounds()
	numcells = reducedsurfacemesh.GetNumberOfCells()
	for p in range(0, numcells):

		ptlist = vtk.vtkIdList()
		normvec = reducedsurfacemesh.GetCellData().GetArray("Normals").GetTuple3(p)

		# If cell normal is in 0,0,1 direction
		if(basalnormal == 'z'):
			if(abs(vtkMath.Dot(normvec, [1,0,0])) < tol and abs(vtkMath.Dot(normvec, [0,1,0])) < tol):
				reducedsurfacemesh.GetCellPoints(p, ptlist)
				count_nodes_at_base = 0;
				for j in range(0, ptlist.GetNumberOfIds()):
					ptid = reducedsurfacemesh.GetPointData().GetArray("vtkIdFilter_Ids").GetValue(ptlist.GetId(j))

					pt = [0,0,0]
					pt = ugrid.GetPoints().GetPoint(ptid)

					if(pt[2] > 0.5*(bds[5] + bds[4])):
						Basalsurfnodes.InsertUniqueId(ptid)

					if(pt[2] > bds[5] - tol):
						count_nodes_at_base = count_nodes_at_base + 1
						
				if(count_nodes_at_base == 3):
					reducedsurfacemesh.DeleteCell(p)

		elif(basalnormal =='x'):

			if(abs(vtkMath.Dot(normvec, [0,0,1])) < tol and abs(vtkMath.Dot(normvec, [0,1,0])) < tol):
				reducedsurfacemesh.GetCellPoints(p, ptlist)
				for j in range(0, ptlist.GetNumberOfIds()):
					ptid = reducedsurfacemesh.GetPointData().GetArray("vtkIdFilter_Ids").GetValue(ptlist.GetId(j))

					pt = [0,0,0]
					pt = ugrid.GetPoints().GetPoint(ptid)


					if(pt[0] < bds[0] + tol):
						Basalsurfnodes.InsertUniqueId(ptid)

					if(pt[0] < bds[0] + tol):
						reducedsurfacemesh.DeleteCell(p)



	reducedsurfacemesh.RemoveDeletedCells();
	#reducedsurfacemesh.Update()

	
	# Split the surfaces to LVendo, RVendo, Epi 
	connectivityfilter = vtk.vtkPolyDataConnectivityFilter()
	connectivityfilter.SetScalarConnectivity(0);
	connectivityfilter.ColorRegionsOn();
	connectivityfilter.SetExtractionModeToAllRegions()
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		connectivityfilter.SetInput(reducedsurfacemesh)
	else:
		connectivityfilter.SetInputData(reducedsurfacemesh)
	connectivityfilter.Update()

	print "Number of connected regions = ", connectivityfilter.GetNumberOfExtractedRegions()
	vtk_py.writeXMLPData(connectivityfilter.GetOutput(), "connectedregion.vtp")


	# Extracting Epi surface
	connectivityfilter.SetExtractionModeToLargestRegion();
	connectivityfilter.ColorRegionsOn();
	connectivityfilter.Update();
	Episurf.DeepCopy(connectivityfilter.GetOutput());


	zrange = []
	for p in range(0, 3):

		connectivityfilter.AddSpecifiedRegion(p);
		connectivityfilter.SetExtractionModeToSpecifiedRegions();
		connectivityfilter.Update();
		connectivityfilter.DeleteSpecifiedRegion(p);

		bds = connectivityfilter.GetOutput().GetBounds()

		if(basalnormal == 'z'):
			zrange.append(abs(bds[5] - bds[4]))
		elif(basalnormal == 'x'):
			zrange.append(abs(bds[1] - bds[0]))

	epiids = zrange.index(max(zrange))
	#RVendoids = zrange.index(min(zrange))
	LVendoids = zrange.index(min(zrange))
	#idlist = set([0,1,2])
	#LVendoids = list(idlist.difference(set([epiids, RVendoids])))[0]

	epiids = 0
	LVendoids = 2
	
	# Extracting RV surface
	[connectivityfilter.AddSpecifiedRegion(p) for p in range(0,connectivityfilter.GetNumberOfExtractedRegions())]
	connectivityfilter.DeleteSpecifiedRegion(epiids);
	connectivityfilter.DeleteSpecifiedRegion(LVendoids);
	connectivityfilter.SetExtractionModeToSpecifiedRegions();
	#connectivityfilter.SetExtractionModeToSpecifiedRegions();
	connectivityfilter.Update();
	RVEndosurf.DeepCopy(connectivityfilter.GetOutput());
	[connectivityfilter.DeleteSpecifiedRegion(p) for p in range(0,connectivityfilter.GetNumberOfExtractedRegions())]
	#connectivityfilter.DeleteSpecifiedRegion(RVendoids);


	outfilename = 'RV_surf.vtp'
	vtk_py.writeXMLPData(RVEndosurf, outfilename)

	# Extracting LV surface
	connectivityfilter.AddSpecifiedRegion(LVendoids);
	connectivityfilter.SetExtractionModeToSpecifiedRegions();
	connectivityfilter.Update();
	LVEndosurf.DeepCopy(connectivityfilter.GetOutput());
	connectivityfilter.DeleteSpecifiedRegion(LVendoids);

	outfilename = 'LV_surf.vtp'
	vtk_py.writeXMLPData(LVEndosurf, outfilename)


	# Get Epi surface points
	for p in range(0, Episurf.GetNumberOfCells()):
		ptlist = vtk.vtkIdList()
		Episurf.GetCellPoints(p, ptlist)
		for j in range(0, ptlist.GetNumberOfIds()):
			ptid = Episurf.GetPointData().GetArray("vtkIdFilter_Ids").GetValue(ptlist.GetId(j))
			Episurfnodes.InsertUniqueId(ptid)


	# Get RV surface points
	for p in range(0, RVEndosurf.GetNumberOfCells()):
		ptlist = vtk.vtkIdList()
		RVEndosurf.GetCellPoints(p, ptlist)
		for j in range(0, ptlist.GetNumberOfIds()):
			ptid = RVEndosurf.GetPointData().GetArray("vtkIdFilter_Ids").GetValue(ptlist.GetId(j))
			RVEndosurfnodes.InsertUniqueId(ptid)

	# Get LV surface points
	for p in range(0, LVEndosurf.GetNumberOfCells()):
		ptlist = vtk.vtkIdList()
		LVEndosurf.GetCellPoints(p, ptlist)
		for j in range(0, ptlist.GetNumberOfIds()):
			ptid = LVEndosurf.GetPointData().GetArray("vtkIdFilter_Ids").GetValue(ptlist.GetId(j))
			LVEndosurfnodes.InsertUniqueId(ptid)

	# Get Epi basal edge points
	Basalepisurfnodes.DeepCopy(Basalsurfnodes)
	Basalepisurfnodes.IntersectWith(Episurfnodes)


	return Basalsurfnodes, Basalepisurfnodes, Episurfnodes, RVEndosurfnodes, LVEndosurfnodes
	


def set_pix_intensity(mesh, vtuoutputdir, vtufilename, txtoutputdir, txtfilename):


	outvtufilename = vtuoutputdir + vtufilename
	outtxtfilename = txtoutputdir + txtfilename

	txtfile = open(outtxtfilename, "w");
	
	cellpixintensity = mesh.GetCellData().GetScalars("closest_pix_intensity")
	print cellpixintensity.GetNumberOfTuples()

	rangeofpixintensity = cellpixintensity.GetValueRange()


	matid = vtk.vtkIntArray()
	matid.SetNumberOfComponents(1)
	matid.SetName("Material Id")

	normalized_pix = vtk.vtkFloatArray()
	normalized_pix.SetNumberOfComponents(1)
	normalized_pix.SetName("Normalized Pixel Intensity")

	for p in range(0, mesh.GetNumberOfCells()):

		pix_intensity = cellpixintensity.GetTuple(p)
		normalized_pix_intensity = (pix_intensity[0] - rangeofpixintensity[0])/(rangeofpixintensity[1] - rangeofpixintensity[0])

		normalized_pix.InsertNextValue(normalized_pix_intensity)

		print >>txtfile, p+1, normalized_pix_intensity

	txtfile.close()
	
	#mesh.GetCellData().SetActiveScalars("Material Id")
	#mesh.GetCellData().SetScalars(matid)

	mesh.GetCellData().SetActiveScalars("Normalized_Pixel_Intensity")
	mesh.GetCellData().SetScalars(normalized_pix)

	writer = vtk.vtkXMLUnstructuredGridWriter()
	writer.SetFileName(outvtufilename)
	writer.SetInput(mesh)
	writer.Write();

def transform_scale_n_write(mesh, outdirectory, vtufilename):

	newfilename = outdirectory + vtufilename[0:len(vtufilename)-4]+"_scaled_rotated.vtu"


	bds = mesh.GetBounds()

	trans = vtk.vtkTransform()
	trans.Translate(bds[5]/10,0,0)
	trans.RotateY(-90)
	trans.Scale(0.1, 0.1, 0.1)


	transfilter = vtk.vtkTransformFilter()
	transfilter.SetTransform(trans)
	transfilter.SetInput(mesh)

	writer = vtk.vtkXMLUnstructuredGridWriter()
	writer.SetFileName(newfilename)
	writer.SetInput(transfilter.GetOutput())
	writer.Write()

def convert_lintets_2_mtv(mesh, outdirectory, name, issavepts, isrotate, basalnormal, tol, epispringcond):


	rotate = isrotate;
	bds = mesh.GetBounds()
	if(basalnormal == 'z'):
		zoffset = bds[5];
		print "z offset = ", zoffset
	elif(basalnormal == 'x'):
		zoffset = bds[1];
		print "x offset = ", zoffset
	
	filename = outdirectory + name + ".grid"
	mtvfile = open(filename, 'w')
	
	print >>mtvfile, " "
	print >>mtvfile, " "
	print >>mtvfile, "Finite element mesh (MLGridFEAdB)"
	print >>mtvfile, " "
	print >>mtvfile, " "
	print >>mtvfile, " "
	print >>mtvfile, "Finite element mesh (GridFE):"
	print >>mtvfile, " "
	print >>mtvfile, "  Number of space dim. =     3  embedded in physical space with dimension 3"
	print >>mtvfile, "  Number of elements   = %5d" %mesh.GetNumberOfCells()
	print >>mtvfile, "  Number of nodes      = %5d" %mesh.GetNumberOfPoints()
	print >>mtvfile, " "
	print >>mtvfile, "  All elements are of the same type : true"
	print >>mtvfile, "  Max number of nodes in an element: 4"
	print >>mtvfile, "  Only one material                : false"
	print >>mtvfile, "  Lattice data                     ? 0"
	print >>mtvfile, " "
	print >>mtvfile, " "
	print >>mtvfile, " "
	print >>mtvfile, " 12 boundary indicators: "
	print >>mtvfile, "   P1"
	print >>mtvfile, "   P2"
	print >>mtvfile, "   P3"
	print >>mtvfile, "   T1"
	print >>mtvfile, "   T2"
	print >>mtvfile, "   u1=0"
	print >>mtvfile, "   u2=0"
	print >>mtvfile, "   u3=0"
	print >>mtvfile, "   u1=u1_0"
	print >>mtvfile, "   u2=u2_0"
	print >>mtvfile, "   u3=u3_0"
	print >>mtvfile, "   free"
	print >>mtvfile, " "
	print >>mtvfile, " "
	print >>mtvfile, " Nodal coordinates and nodal boundary indicators,"
	print >>mtvfile, " the columns contain:"
	print >>mtvfile, "  - node number"
	print >>mtvfile, "  - coordinates"
	print >>mtvfile, "  - no of boundary indicators that are set (ON)"
	print >>mtvfile, "  - the boundary indicators that are set (ON) if any."
	print >>mtvfile, "#"


	# Get Surface nodes for pressure constraint
	Basalsurfnodes, Basalepisurfnodes, Episurfnodes, RVEndosurfnodes, LVEndosurfnodes = extract_biven_surf_nodes(mesh, basalnormal, tol)

	if(issavepts == 1):
		ptfilename = outdirectory + name + '_Basalsurfnodes.vtp'
		savepoints(mesh, Basalsurfnodes, ptfilename)

		ptfilename = outdirectory + name + '_Basalepisurfnodes.vtp'
		savepoints(mesh, Basalepisurfnodes, ptfilename)

		ptfilename = outdirectory + name + '_Episurfnodes.vtp'
		savepoints(mesh, Episurfnodes, ptfilename)

		ptfilename = outdirectory + name + '_RVEndosurfnodes.vtp'
		savepoints(mesh, RVEndosurfnodes, ptfilename)

		ptfilename = outdirectory + name + '_LVEndosurfnodes.vtp'
		savepoints(mesh, LVEndosurfnodes, ptfilename)



	
	num_bindicator = zeros(mesh.GetNumberOfPoints());
	bindicator = ['']*mesh.GetNumberOfPoints();
	for k in range(0,Episurfnodes.GetNumberOfIds()):
		nodeid = Episurfnodes.GetId(k)
		num_bindicator[nodeid] = num_bindicator[nodeid] + 1
		if(epispringcond==1):
			bindicator[nodeid] = bindicator[nodeid] + ' 4'
		elif(epispringcond==2):
			bindicator[nodeid] = bindicator[nodeid] + ' 4'
		else:
			bindicator[nodeid] = bindicator[nodeid] + ' 3'

	for k in range(0,LVEndosurfnodes.GetNumberOfIds()):
		nodeid = LVEndosurfnodes.GetId(k)
		num_bindicator[nodeid] = num_bindicator[nodeid] + 1
		bindicator[nodeid] = bindicator[nodeid] + ' 1'

	for k in range(0,RVEndosurfnodes.GetNumberOfIds()):
		nodeid = RVEndosurfnodes.GetId(k)
		num_bindicator[nodeid] = num_bindicator[nodeid] + 1
		bindicator[nodeid] = bindicator[nodeid] + ' 2'

	for k in range(0,Basalepisurfnodes.GetNumberOfIds()):
		nodeid = Basalepisurfnodes.GetId(k)
		if(epispringcond==1):
			num_bindicator[nodeid] = num_bindicator[nodeid] + 0
		elif(epispringcond==2):
			num_bindicator[nodeid] = num_bindicator[nodeid] + 2
			if(rotate == 0):
				bindicator[nodeid] = bindicator[nodeid] + ' 4 8'
			else:
				bindicator[nodeid] = bindicator[nodeid] + ' 4 6'
		else:
			num_bindicator[nodeid] = num_bindicator[nodeid] + 3
			bindicator[nodeid] = bindicator[nodeid] + ' 6 7 8'

	for k in range(0,Basalsurfnodes.GetNumberOfIds()):
		nodeid = Basalsurfnodes.GetId(k)
		if(rotate == 0):
			if(Basalepisurfnodes.IsId(nodeid) == -1):
				if(epispringcond==1):
					num_bindicator[nodeid] = num_bindicator[nodeid] + 0
				elif(epispringcond==2):
					num_bindicator[nodeid] = num_bindicator[nodeid] + 1
					bindicator[nodeid] = bindicator[nodeid] + ' 8'
				else:
					num_bindicator[nodeid] = num_bindicator[nodeid] + 1
					bindicator[nodeid] = bindicator[nodeid] + ' 8'
		else:
			if(Basalepisurfnodes.IsId(nodeid) == -1):
				if(epispringcond==1):
					num_bindicator[nodeid] = num_bindicator[nodeid] + 0
				elif(epispringcond==2):
					num_bindicator[nodeid] = num_bindicator[nodeid] + 1
					bindicator[nodeid] = bindicator[nodeid] + ' 6'
				else:
					num_bindicator[nodeid] = num_bindicator[nodeid] + 1
					bindicator[nodeid] = bindicator[nodeid] + ' 6'
	
	#print bindicator
	pt = [0,0,0];
	rotatedpt = [0,0,0];
	#num_bindicator = 0;
	#bindicator = " ";
	for p in range(0, mesh.GetNumberOfPoints()):
		if(rotate == 0):
			mesh.GetPoints().GetPoint(p,pt);
	     		print >>mtvfile, "%6d  ( %11.5e, %11.5e, %11.5e)  [%d] %s " %((p+1), pt[0], pt[1], pt[2], num_bindicator[p], bindicator[p] )
		else:
			mesh.GetPoints().GetPoint(p,pt);
			rotatedpt[0] = (-pt[2]+zoffset)/10.0;
			rotatedpt[1] = (pt[1])/10.0;
			rotatedpt[2] = (pt[0])/10.0;
			
	     		print >>mtvfile, "%6d  ( %11.5e, %11.5e, %11.5e)  [%d] %s " %((p+1), rotatedpt[0], rotatedpt[1], rotatedpt[2], num_bindicator[p], bindicator[p] )
	
	print >>mtvfile, " "
	print >>mtvfile, "  Element types and connectivity"
	print >>mtvfile, "  the columns contain:"
	print >>mtvfile, "   - element number"
	print >>mtvfile, "   - element type"
	print >>mtvfile, "   - material number"
	print >>mtvfile, "   - the global node numbers of the nodes in the element."
	print >>mtvfile, "#"


	materialid = mesh.GetCellData().GetArray('Material Id')
	ptid = vtkIdList()
	ids = [0,0,0,0,0,0,0,0,0,0];
	
	for p in range(0, mesh.GetNumberOfCells()):
		mesh.GetCellPoints(p, ptid);
		for j in range(0, 4):
			ids[j] = int(ptid.GetId(j));
			try:
				matid = materialid.GetValue(p)

			except AttributeError:
				matid = 1

	     	print >>mtvfile, "%5d  ElmT4n3D  %1d       %5d %5d %5d %5d " %((p+1), matid, ids[0]+1, ids[1]+1, ids[2]+1, ids[3]+1)

		
	mtvfile.close()



tol = 1e-1
parser = argparse.ArgumentParser()
parser.add_argument('--vtk_folder', type=str, required=True)
parser.add_argument('--vtk_filename', type=str, required=True)
parser.add_argument('--mtv_grid_directory', type=str, required=True)
parser.add_argument('--mtv_basename', type=str, required=True)
parser.add_argument('--isepispring', type=int, required=True)
args = parser.parse_args()


print "************* Entering vtkmtvtranslator_lintets_v1.py *****************"

if(args.isepispring):
	print "Enforce spring B.C on epi"


if (args.vtk_filename[len(args.vtk_filename)-3:len(args.vtk_filename)] == 'vtu'):
	mesh = vtk_py.readXMLUGrid(os.path.join(args.vtk_folder, args.vtk_filename))
else:
	mesh = vtk_py.readUGrid(os.path.join(args.vtk_folder, args.vtk_filename))

mesh = extract_tetra(mesh)

convert_lintets_2_mtv(mesh, args.mtv_grid_directory, args.mtv_basename, 1, 1, "z", tol, args.isepispring)

print "************* Leaving vtkmtvtranslator_lintets_v1.py *****************"
