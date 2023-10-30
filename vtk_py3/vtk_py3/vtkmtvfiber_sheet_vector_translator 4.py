########################################################################

import argparse
import glob
from numpy import *
from sets import Set
from vtk import *
import os

import vtk_py as vtk_py

########################################################################



parser = argparse.ArgumentParser()
parser.add_argument('--vtk_folder', type=str, required=True)
parser.add_argument('--vtk_filename', type=str, required=True)
parser.add_argument('--fiberdata_directory', type=str, required=True)
parser.add_argument('--mtv_grid_directory', type=str, required=True)
parser.add_argument('--mtv_basename', type=str, required=True)
args = parser.parse_args()

print "************* Entering vtkmtvfiber_sheet_vector_translator.py *****************"


points = []
fvectors = []
svectors = []

srcfilename = os.path.join(args.vtk_folder, args.vtk_filename)
directory = args.fiberdata_directory
outdirectory = args.mtv_grid_directory
outfilename = args.mtv_basename

reader = vtkUnstructuredGridReader();
reader.SetFileName(srcfilename);
reader.Update();

mesh = vtkUnstructuredGrid();
mesh.DeepCopy(reader.GetOutput())

rotate = 1;
bds = mesh.GetBounds()
zoffset = bds[5];


rotate = 1;
# Import vtk file
# Read in LV fiber data
LVfiberreader = vtkXMLPolyDataReader();
LVfiberreader.SetFileName(directory+"LVfiberdata.vtp");
LVfiberreader.Update();

LVfiberdata = vtkPolyData();
LVfiberdata.DeepCopy(LVfiberreader.GetOutput())

fvec = [0,0,0]
svec = [0,0,0]
pts = [0,0,0]
rotatedpt = [0,0,0]

for p in range(0, LVfiberdata.GetPoints().GetNumberOfPoints()):
	pt= LVfiberdata.GetPoints().GetPoint(p)
	if(rotate == 0):	
		points.append(pt)
	else:
		rotatedpt[0] = float(-pt[2]+zoffset)/10.0;
		rotatedpt[1] = float(pt[1])/10.0;
		rotatedpt[2] = float(pt[0])/10.0;
		points.append([rotatedpt[0], rotatedpt[1], rotatedpt[2]])


	LVfiberdata.GetPointData().SetActiveVectors("f vectors")
	temp = LVfiberdata.GetPointData().GetVectors().GetTuple3(p);
	fvec[0] = temp[0]
	fvec[1] = temp[1]
	fvec[2] = temp[2]
	vtkMath.Normalize(fvec)
	if(rotate == 0):
		fvectors.append([fvec[0], fvec[1], fvec[2]])
	else:
		fvectors.append([-1.0*fvec[2], fvec[1], fvec[0]])


	LVfiberdata.GetPointData().SetActiveVectors("n vectors")
	temp = LVfiberdata.GetPointData().GetVectors().GetTuple3(p);
	svec[0] = temp[0]
	svec[1] = temp[1]
	svec[2] = temp[2]
	vtkMath.Normalize(svec)
	if(rotate == 0):
		svectors.append([svec[0], svec[1], svec[2]])
	else:
		svectors.append([-1.0*svec[2], svec[1], svec[0]])
		

# Read in RV fiber data
RVfiberreader = vtkXMLPolyDataReader();
RVfiberreader.SetFileName(directory+"RVfiberdata.vtp");
RVfiberreader.Update();

RVfiberdata = vtkPolyData();
RVfiberdata.DeepCopy(RVfiberreader.GetOutput())

fvec = [0,0,0]
svec = [0,0,0]
pts = [0,0,0]

rotatedpt = [0,0,0]

for p in range(0, RVfiberdata.GetPoints().GetNumberOfPoints()):
	pt= RVfiberdata.GetPoints().GetPoint(p)
	if(rotate == 0):	
		points.append(pt)
	else:
		rotatedpt[0] = float(-pt[2]+zoffset)/10.0;
		rotatedpt[1] = float(pt[1])/10.0;
		rotatedpt[2] = float(pt[0])/10.0;
		points.append([rotatedpt[0], rotatedpt[1], rotatedpt[2]])



	RVfiberdata.GetPointData().SetActiveVectors("f vectors")
	temp = RVfiberdata.GetPointData().GetVectors().GetTuple3(p);
	fvec[0] = temp[0]
	fvec[1] = temp[1]
	fvec[2] = temp[2]
	vtkMath.Normalize(fvec)
	if(rotate == 0):
		fvectors.append([fvec[0], fvec[1], fvec[2]])
	else:
		fvectors.append([-1.0*fvec[2], fvec[1], fvec[0]])

	RVfiberdata.GetPointData().SetActiveVectors("n vectors")
	temp = RVfiberdata.GetPointData().GetVectors().GetTuple3(p);
	svec[0] = temp[0]
	svec[1] = temp[1]
	svec[2] = temp[2]
	vtkMath.Normalize(svec)
	if(rotate == 0):
		svectors.append([svec[0], svec[1], svec[2]])
	else:
		svectors.append([-1.0*svec[2], svec[1], svec[0]])



# Read in Septum fiber data
Septumfiberreader = vtkXMLPolyDataReader();
Septumfiberreader.SetFileName(directory+"Septumfiberdata.vtp");
Septumfiberreader.Update();

Septumfiberdata = vtkPolyData();
Septumfiberdata.DeepCopy(Septumfiberreader.GetOutput())

fvec = [0,0,0]
svec = [0,0,0]
pts = [0,0,0]
rotatedpt = [0,0,0]


for p in range(0, Septumfiberdata.GetPoints().GetNumberOfPoints()):
	pt= Septumfiberdata.GetPoints().GetPoint(p)
	if(rotate == 0):	
		points.append(pt)
	else:
		rotatedpt[0] = float(-pt[2]+zoffset)/10.0;
		rotatedpt[1] = float(pt[1])/10.0;
		rotatedpt[2] = float(pt[0])/10.0;
		points.append([rotatedpt[0], rotatedpt[1], rotatedpt[2]])


	Septumfiberdata.GetPointData().SetActiveVectors("f vectors")
	temp = Septumfiberdata.GetPointData().GetVectors().GetTuple3(p);
	fvec[0] = temp[0]
	fvec[1] = temp[1]
	fvec[2] = temp[2]
	vtkMath.Normalize(fvec)
	if(rotate == 0):
		fvectors.append([fvec[0], fvec[1], fvec[2]])
	else:
		fvectors.append([-1.0*fvec[2], fvec[1], fvec[0]])


	Septumfiberdata.GetPointData().SetActiveVectors("n vectors")
	temp = Septumfiberdata.GetPointData().GetVectors().GetTuple3(p);
	svec[0] = temp[0]
	svec[1] = temp[1]
	svec[2] = temp[2]
	vtkMath.Normalize(svec)
	if(rotate == 0):
		svectors.append([svec[0], svec[1], svec[2]])
	else:
		svectors.append([-1.0*svec[2], svec[1], svec[0]])
		

mtvfiberfilename = outdirectory+outfilename+"_fiber_rotated.axis"
mtvfiberfile = open(mtvfiberfilename, 'w')
print >>mtvfiberfile, len(fvectors), 10

mtvsheetfilename = outdirectory+outfilename+"_sheet_rotated.axis"
mtvsheetfile = open(mtvsheetfilename, 'w')
print >>mtvsheetfile, len(svectors)

pdata = vtk.vtkPolyData()
ppoints = vtk.vtkPoints()
cellarray = vtk.vtkCellArray()

fvector = vtk.vtkFloatArray()
fvector.SetNumberOfComponents(3)
fvector.SetName("f vectors")
svector = vtk.vtkFloatArray()
svector.SetNumberOfComponents(3)
svector.SetName("s vectors")
nvector = vtk.vtkFloatArray()
nvector.SetNumberOfComponents(3)
nvector.SetName("n vectors")

for p in range(0,len(fvectors)):

	print >>mtvfiberfile, points[p][0], points[p][1], points[p][2], fvectors[p][0], fvectors[p][1], fvectors[p][2]
	print >>mtvsheetfile, points[p][0], points[p][1], points[p][2], svectors[p][0], svectors[p][1], svectors[p][2]

	ppoints.InsertNextPoint(points[p][0],points[p][1], points[p][2])

	vert = vtk.vtkVertex()
	vert.GetPointIds().SetId(0,p)
	cellarray.InsertNextCell(vert)

	fvector.InsertNextTuple3(fvectors[p][0], fvectors[p][1], fvectors[p][2])
	svector.InsertNextTuple3(svectors[p][0], svectors[p][1], svectors[p][2])

	f = array([fvectors[p][0], fvectors[p][1], fvectors[p][2]])  
	s = array([svectors[p][0], svectors[p][1], svectors[p][2]])  
	n = cross(f,s)
	
	nvector.InsertNextTuple3(n[0], n[1], n[2])

pdata.SetPoints(ppoints)
pdata.SetVerts(cellarray)
pdata.GetPointData().SetActiveVectors("f vectors")
pdata.GetPointData().SetVectors(fvector)
pdata.GetPointData().SetActiveVectors("s vectors")
pdata.GetPointData().SetVectors(svector)
pdata.GetPointData().SetActiveVectors("n vectors")
pdata.GetPointData().SetVectors(nvector)


newfilename = outdirectory + outfilename + "_scaled_rotated.vtp"

writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(newfilename)
if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
	writer.SetInput(pdata)
else:
	writer.SetInputData(pdata)
writer.Write()

mtvfiberfile.close()
mtvsheetfile.close()


print "************* Leaving vtkmtvfiber_sheet_vector_translator.py *****************"
