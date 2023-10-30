########################################################################

import numpy
import vtk

########################################################################

def readMSHGrid(mshfilename):

	fdata = open(mshfilename, "r")
	lines = fdata.readlines()
	startelem = 0
	startnode = 0
	cnt = 0

	cellArray = vtk.vtkCellArray()
	points = vtk.vtkPoints()
	for line in lines:
		if("$Elements" in line):
			startelem = 1
			continue

		if("$EndElements" in line):
			startelem = 0
			continue

		if("$Nodes" in line):
			startnode = 1
			continue

		if("$EndNodes" in line):
			startnode = 0
			continue

		if(startnode == 1):
			try:
				x = float(line.strip().split()[1])
				y = float(line.strip().split()[2])
				z = float(line.strip().split()[3])

				points.InsertNextPoint(x, y, z)

			except IndexError:
				continue


		if(startelem == 1):
			try:
				elmtype = int(line.strip().split()[1])
				if(elmtype == 4):
					tetra = vtk.vtkTetra()
					tetra.GetPointIds().SetId(0, int(line.strip().split()[5])-1)
					tetra.GetPointIds().SetId(1, int(line.strip().split()[6])-1)
					tetra.GetPointIds().SetId(2, int(line.strip().split()[7])-1)
					tetra.GetPointIds().SetId(3, int(line.strip().split()[8])-1)
					cellArray.InsertNextCell(tetra)
			except IndexError:
				continue

	ugrid = vtk.vtkUnstructuredGrid()
	ugrid.SetPoints(points)
	ugrid.SetCells(10, cellArray)

	return ugrid



