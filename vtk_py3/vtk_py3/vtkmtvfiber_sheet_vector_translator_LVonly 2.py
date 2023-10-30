########################################################################

import argparse
import glob
from numpy import *
from sets import Set
from vtk import *
import os

import vtk_py as vtk_py

########################################################################

def transform_scale_n_write(mesh, outdirectory):

	newfilename = outdirectory + "fiberdata" +"_scaled_rotated.vtp"

	bds = mesh.GetBounds()

	trans = vtk.vtkTransform()
	trans.Translate(bds[5]/10,0,0)
	trans.RotateY(-90)
	trans.Scale(0.1, 0.1, 0.1)

	transfilter = vtk.vtkTransformPolyDataFilter()
	transfilter.SetTransform(trans)
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		transfilter.SetInput(mesh)
	else:
		transfilter.SetInputData(mesh)
	transfilter.Update()

	writer = vtk.vtkXMLPolyDataWriter()
	writer.SetFileName(newfilename)
	if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
		writer.SetInput(transfilter.GetOutput())
	else:
		writer.SetInputData(transfilter.GetOutput())
	writer.Write()

parser = argparse.ArgumentParser()
parser.add_argument('--vtk_folder', type=str, required=True)
parser.add_argument('--vtk_filename', type=str, required=True)
parser.add_argument('--fiberdata_directory', type=str, required=True)
parser.add_argument('--mtv_grid_directory', type=str, required=True)
parser.add_argument('--mtv_basename', type=str, required=True)
parser.add_argument('--scale', type=float, required=True)
args = parser.parse_args()

print "************* Entering vtkmtvfiber_sheet_vector_translator_LVonly.py *****************"


points = []
fvectors = []
svectors = []

srcfilename = os.path.join(args.vtk_folder, args.vtk_filename)
directory = args.fiberdata_directory
outdirectory = args.mtv_grid_directory
outfilename = args.mtv_basename
scale = args.scale

#srcfilename = "/home/likchuan/Research/C++/heartmech_bitbucket/pulse_cmake/pulse-cmake-growth/cases/singleLV/mesh/LV_fine.vtu"
#directory = "/home/likchuan/Research/C++/heartmech_bitbucket/pulse_cmake/pulse-cmake-growth/cases/singleLV/mesh/"
#outdirectory = "/home/likchuan/Research/C++/heartmech_bitbucket/pulse_cmake/pulse-cmake-growth/cases/singleLV/mesh/"
#outfilename = "LV_fiber_fine"

reader = vtkXMLUnstructuredGridReader();
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
LVfiberreader.SetFileName(directory+"fiberdirections.vtp");
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
		rotatedpt[0] = float(-pt[2]+zoffset)/scale;
		rotatedpt[1] = float(pt[1])/scale;
		rotatedpt[2] = float(pt[0])/scale;
		#print rotatedpt
		points.append([rotatedpt[0], rotatedpt[1], rotatedpt[2]])


	LVfiberdata.GetPointData().SetActiveVectors("fiber vectors")
	temp = LVfiberdata.GetPointData().GetVectors().GetTuple3(p);
	fvec[0] = temp[0]
	fvec[1] = temp[1]
	fvec[2] = temp[2]
	vtkMath.Normalize(fvec)
	if(rotate == 0):
		fvectors.append([fvec[0], fvec[1], fvec[2]])
	else:
		fvectors.append([-1.0*fvec[2], fvec[1], fvec[0]])


	LVfiberdata.GetPointData().SetActiveVectors("sheet vectors")
	temp = LVfiberdata.GetPointData().GetVectors().GetTuple3(p);
	svec[0] = temp[0]
	svec[1] = temp[1]
	svec[2] = temp[2]
	vtkMath.Normalize(svec)
	if(rotate == 0):
		svectors.append([svec[0], svec[1], svec[2]])
	else:
		svectors.append([-1.0*svec[2], svec[1], svec[0]])
		

# Appending PolyData
appendeddata = vtk.vtkAppendPolyData()
if(vtk.vtkVersion.GetVTKMajorVersion() < 6):
	appendeddata.AddInput(LVfiberreader.GetOutput())
else:
	appendeddata.AddInputData(LVfiberreader.GetOutput())
appendeddata.Update()

#print appendeddata.GetOutput()

#transform_scale_n_write(appendeddata.GetOutput(), outdirectory)


mtvfiberfilename = outdirectory+outfilename+"_fiber_rotated.axis"
mtvfiberfile = open(mtvfiberfilename, 'w')
print >>mtvfiberfile, len(fvectors), 10

mtvsheetfilename = outdirectory+outfilename+"_sheet_rotated.axis"
mtvsheetfile = open(mtvsheetfilename, 'w')
print >>mtvsheetfile, len(svectors)

for p in range(0,len(fvectors)):
	print >>mtvfiberfile, points[p][0], points[p][1], points[p][2], fvectors[p][0], fvectors[p][1], fvectors[p][2]
	print >>mtvsheetfile, points[p][0], points[p][1], points[p][2], svectors[p][0], svectors[p][1], svectors[p][2]


mtvfiberfile.close()
mtvsheetfile.close()


print "************* Leaving vtkmtvfiber_sheet_vector_translator.py *****************"
